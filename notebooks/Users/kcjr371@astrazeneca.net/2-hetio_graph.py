# Databricks notebook source
# MAGIC %md
# MAGIC 
# MAGIC # Creating graph and calculating metapaths
# MAGIC 
# MAGIC __Purpose__: To calculate the metapaths for each compound/disease pair in the het.io graph
# MAGIC 
# MAGIC __Input__: Two csv files representing nodes and edges (/zenegraph/bikg-sources/hetio_nodes.csv and /zenegraph/bikg-sources/hetio_edges.csv)
# MAGIC 
# MAGIC __Output__: ? not clear yet...

# COMMAND ----------

from pyspark.sql import *
from pyspark.sql.functions import *
from pyspark.sql.types import *
from functools import reduce
import pandas as pd
import itertools
import copy
#from pyspark.sql.functions import array, udf
#from pyspark.sql.types import ArrayType, StringType, IntegerType

# COMMAND ----------

# enabling access to the blob storage
storage_account_name = "zenegraph"
storage_account_access_key = dbutils.secrets.get(scope = "zenegraph", key = "zenegraph-storage")

spark.conf.set(
  "fs.azure.account.key."+storage_account_name+".blob.core.windows.net",
  storage_account_access_key)
# needed to use pandas
spark.conf.set("spark.sql.execution.arrow.enabled", "true")

# COMMAND ----------

# MAGIC %md
# MAGIC ## Loading the data

# COMMAND ----------

# function for loading csv files from blob storage
def read_from_azure_blob(filename):
  read_location = "wasbs://bikg-sources@zenegraph.blob.core.windows.net/"+filename+".csv"
  df = spark.read.format("csv").option("header", "true").load(read_location)
  return df
 
# load node csv + rename columns
nodes = read_from_azure_blob("hetio_nodes")
nodes = nodes.selectExpr('id as id', 'identifiers as identifier', 'labels as description', 'type as entity_type')

# load edge csv + rename columns
edges_original = read_from_azure_blob("hetio_edges")
edges_original = edges_original.selectExpr('src as src', 'dst as dst', 'type as edge_type')
edges_original.count()

# COMMAND ----------

# MAGIC %md
# MAGIC # Het.io graph

# COMMAND ----------

# MAGIC %md
# MAGIC ### Step 1: Create bi-directional edges
# MAGIC 
# MAGIC In the original het.io graph, each edge represents a bidrectional relationship. GraphFrames does not accept undirected relationships so it thinks that each edge I have put in the edges dataframe, is a directed one. I will add the opposite edge to the dataframe so that I can then calculate paths. Or will this mess it up?!?!

# COMMAND ----------

edges_temp = edges_original
edges_extra = edges_temp.selectExpr('dst as src', 'src as dst', 'edge_type as edge_type')

# append edges_extra to edges
import functools
def unionAll(dfs):
    return functools.reduce(lambda df1,df2: df1.union(df2.select(df1.columns)), dfs) 

edges = unionAll([edges_original, edges_extra])
edges.count()

# COMMAND ----------

# MAGIC %md
# MAGIC 
# MAGIC ### Step 2: Create the graph
# MAGIC 
# MAGIC Using the original nodes, and the updated edges (to imitate an undirected graph)

# COMMAND ----------

from graphframes import *

g = GraphFrame(nodes, edges)
print(g)

# COMMAND ----------

# MAGIC %md
# MAGIC 
# MAGIC Calculate all inDegrees - This should correspond to the number of degrees in the original het.io graph (recall that I have added extra degrees to be able to calculate the metapaths)

# COMMAND ----------

node_degrees = g.inDegrees
node_degrees.show()

# COMMAND ----------

# MAGIC %md
# MAGIC ## Find a metapath
# MAGIC 
# MAGIC 1. Select a metapath: Compound-[binds]-(Gene)-[associated]-(Disease)
# MAGIC 2. Select a compound/disease pair: bupropion (DB01156)/nicotine dependence (DOID:0050742)
# MAGIC 2. Calculate all paths through that metapath for the specific pair
# MAGIC 3. Calculate DWPC

# COMMAND ----------

def metapaths_of_length_2(origin_node, r1_type, middle_node, r2_type, destination_node):
  motifs = g.find("(a)-[r1]->(b); (b)-[r2]->(c)") 
  metapaths_df = motifs.filter("a.identifier ='"+origin_node+"' and r1.edge_type = '"+r1_type+"' and b.entity_type='"+middle_node+"' and r2.edge_type='"+r2_type+"' and c.identifier='"+destination_node+"'")
  return metapaths_df

# COMMAND ----------

mp_2 = metapaths_of_length_2('DB01156', 'binds', 'Gene', 'associates', 'DOID:0050742')
display(mp_2)

# COMMAND ----------

# MAGIC %md
# MAGIC I now need to extract the node ids from the dataframe mp_2

# COMMAND ----------



# COMMAND ----------

temp = mp_2.select('a').collect()

# COMMAND ----------

temp.id

# COMMAND ----------

# MAGIC %md
# MAGIC Calculating the number of edge types

# COMMAND ----------

# Find chains of 4 vertices.
chain4 = g.find("(a)-[ab]->(b); (b)-[bc]->(c); (c)-[cd]->(d)")

# Query on sequence, with state (cnt)
#  (a) Define method for updating state given the next element of the motif.
def path_degree_products(node):
  g.degrees.filter('node')
  return when(relationship == "friend", cnt + 1).otherwise(cnt)

#  (b) Use sequence operation to apply method to sequence of elements in motif.
#   In this case, the elements are the 3 edges.
edges = ["ab", "bc", "cd"]
numFriends = reduce(cumFriends, edges, lit(0))
    
chainWith2Friends2 = chain4.withColumn("num_friends", numFriends).where(numFriends >= 2)
display(chainWith2Friends2)