# Databricks notebook source
# MAGIC %md
# MAGIC 
# MAGIC # Creating graph and calculating metapaths
# MAGIC 
# MAGIC __Purpose__: To calculate the metapaths for each compound/disease pair in the het.io graph
# MAGIC 
# MAGIC __Input__: Two csv files representing nodes and edges (/zenegraph/bikg-sources/hetio_nodes.csv and /zenegraph/bikg-sources/hetio_edges.csv)
# MAGIC 
# MAGIC __Output__: Dataframe with three columns: start_node, end_node and DWPC

# COMMAND ----------

from pyspark.sql import *
from pyspark.sql.functions import *
from pyspark.sql.types import *
from functools import reduce
import pandas as pd
import itertools
import copy

from pyspark.sql.functions import col

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

# The number of unique compounds in the graph
nodes.filter("entity_type = 'Compound'").select(col('identifier')).count()

# COMMAND ----------

# The number of unique diseases in the graph
nodes.filter("entity_type = 'Disease'").select(col('identifier')).count()

# COMMAND ----------

1552*137

# COMMAND ----------

# MAGIC %md
# MAGIC # Het.io graph

# COMMAND ----------

# MAGIC %md
# MAGIC ## Step 1: Create bi-directional edges
# MAGIC 
# MAGIC In the original het.io graph, each edge represents a bidrectional relationship. GraphFrames does not accept undirected relationships so it thinks that each edge I have put in the edges dataframe, is a directed one. I will add the opposite edge to the dataframe so that I can then calculate paths. 

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
# MAGIC ## Step 2: Create the graph
# MAGIC 
# MAGIC Using the original nodes, and the updated edges (to imitate an undirected graph)

# COMMAND ----------

from graphframes import *

g = GraphFrame(nodes, edges)
print(g)

# COMMAND ----------

# MAGIC %md
# MAGIC The following extracts all the metapaths of length 3 for a specified compound/disease pair

# COMMAND ----------

g_original = GraphFrame(nodes, edges_original)

# metapaths of length 2
metapath = g_original.find("(c)-[r1]->(n1); (n1)-[r2]->(d)")

# filter for the specific compound and disease
metapath_filter = metapath.filter("c.identifier ='{}' and d.identifier='{}'".format('DB01156', 'DOID:0050742'))

# extract relevant information from each column (e.g., node and edge types)
metapath_df = metapath_filter.withColumn('r1', col('r1').getField('edge_type')).withColumn('n1', col('n1').getField('entity_type')).withColumn('r2', col('r2').getField('edge_type'))

# aggregate the dataframe by metapath and count how many paths there are per metapath
metapath_df.groupby(['c','r1', 'n1', 'r2','d']).agg({'d':'count'}).withColumnRenamed('count(d)', 'number_of_paths').sort(desc('count(d)')).show()

# COMMAND ----------

metapath_filter.show()

# COMMAND ----------

# MAGIC %md
# MAGIC ## Step 3: Find a metapath
# MAGIC 
# MAGIC 1. Select a metapath: e.g., Compound-[binds]-(Gene)-[associated]-(Disease)
# MAGIC 
# MAGIC 2. Select a compound/disease pair: e.g., bupropion (DB01156)/nicotine dependence (DOID:0050742)
# MAGIC 
# MAGIC 2. Extract all paths along the specified metapath, for the specified compound/disease pair
# MAGIC 
# MAGIC 3. Calculate DWPC
# MAGIC 
# MAGIC   3.1. Select the columns with '\_node' as a suffix
# MAGIC   
# MAGIC   3.2. Extract 'id' from each column and replace current content with only this 'id'
# MAGIC   
# MAGIC   3.3. Get the node degrees (based on the node 'id') and multiply these degrees by 1**-0.4 (as per the rephetio analysis)
# MAGIC   
# MAGIC   3.4. Create a new column as the sum of all columns, call this column 'pdp'
# MAGIC   
# MAGIC   3.5. Collapse or groupby 'start_node' and 'end_node' and sum (pdp) -> this is the DWPC for the specific metapath and start/end nodes

# COMMAND ----------

# MAGIC %md
# MAGIC Helper functions below
# MAGIC 1. metapaths_of_length_n (for n=2,3,4) takes as input the elements along the path and returns a dataframe where each row represents a path along that metapath
# MAGIC 2. column_add takes as input two columns and returns their sum

# COMMAND ----------

# helper functions/dataframes

# this function calculates the paths along the specific metapath and returns them in a dataframe, each column represents an element along the path - node or edge
def metapaths_of_length_2(compound, r1, n1, r2, disease):  
  motifs = g.find("(compound)-[r1]->(n1); (n1)-[r2]->(disease)") 
  metapaths_df = motifs.filter("compound.identifier ='{}' and r1.edge_type = '{}' and n1.entity_type='{}' and r2.edge_type='{}' and disease.identifier='{}'".format(compound, r1, n1, r2, disease))
  return metapaths_df

# metapaths of length 3
def metapaths_of_length_3(compound, r1, n1, r2, n2, r3, disease):  
  motifs = g.find("(compound)-[r1]->(n1); (n1)-[r2]->(n2); (n2)-[r3]->(disease)") 
  metapaths_df = motifs.filter("compound.identifier ='{}' and r1.edge_type = '{}' and n1.entity_type='{}' and r2.edge_type='{}' and n2.entity_type='{}' and r3.edge_type='{}' and disease.identifier='{}'".format(compound, r1, n1, r2, n2, r3, disease))
  return metapaths_df

# metapaths of length 4
def metapaths_of_length_4(compound, r1, n1, r2, n2, r3, n3, r4, disease):  
  motifs = g.find("(compound)-[r1]->(n1); (n1)-[r2]->(n2); (n2)-[r3]->(n3); (n3)-[r4]-(disease)") 
  metapaths_df = motifs.filter("compound.identifier ='{}' and r1.edge_type = '{}' and n1.entity_type='{}' and r2.edge_type='{}' and n2.entity_type='{}' and r3.edge_type='{}' and n3.entity_type='{}' and r3.edge_type='{}' and disease.identifier='{}'".format(compound, r1, n1, r2, n2, r3, n3, r4, disease))
  return metapaths_df

# this function is for when I want to add all the columns for the pdp
def column_add(a,b):
     return  a.__add__(b) 
  
# dataframe of in-degrees - this should correspond to the original het.io degrees for each node
node_degrees = g.inDegrees

# COMMAND ----------

# MAGIC %md
# MAGIC DWPC function, returns dataframe with DWPC for specified compound/disease pair and specific metapath

# COMMAND ----------

# this function takes as input the specified metapath (each element of the metapath is an input) along with a dataframe of the in-degrees of each node

def dwpc_2(compound, r1, n1, r2, disease, degrees_df):

  # step 0: calculate metapath dataframe  
  df = metapaths_of_length_2(compound, r1, n1, r2, disease)
  
  # step 1: keep node columns
  condition = lambda col: '_node' not in col
  df = df.drop(*filter(condition, df.columns))
  
  # step 2: only keep 'id' value in each column
  for c in df.columns:
    df = df.withColumn(c, col(c).getField('id'))
  
  # step 3: Get the degrees and multiply by 1**-0.4. Probably a much better way of doing this... 
  columns = df.columns
  for c in columns:
    temp1 = df.withColumnRenamed(c, 'id')
    temp2 = temp1.join(node_degrees, 'id')
    temp3 = temp2.withColumnRenamed('inDegree', c).drop('id').withColumn(c, col(c)**-0.4) 
    df = temp3

  # step 4: sum all the columns to obtain the pdp
  pdp_df = df.withColumn('pdp', reduce(column_add, (df[col] for col in df.columns)))
  
  # step 5: replace start_node and end_node columns with their appropriate ids
  pdp_df = pdp_df.withColumn('compound', lit(start_node)).withColumn('disease', lit(end_node))
  
  # step 6: group by start_node and end_node to obtain DWPC
  dwpc_df = pdp_df.groupby(['compound', 'disease']).agg({'pdp':'sum'}).withColumnRenamed('sum(pdp)', 'dwpc')
  
  return dwpc_df

# COMMAND ----------

dwpc_2('DB01156', 'binds', 'Gene', 'associates', 'DOID:0050742', node_degrees).show()

# COMMAND ----------



# COMMAND ----------

