# Databricks notebook source
# MAGIC %md
# MAGIC # Converting het.io json files into GraphFrames appropriate format
# MAGIC 
# MAGIC __Purpose__: This notebook takes two json files of node and edge information and converts them to the appropriate format for GraphFrames. 
# MAGIC 
# MAGIC __Input__: Two json files with node and edge information from the het.io graph. These files are contained in the blob storage as zenegraph/bikg-sources/hetionet_edges.json and zenegraph/bikg-sources/hetionet_nodes.json
# MAGIC 
# MAGIC __Output__: Two csv files stored in zenegraph/bikg-sources, from two dataframes:
# MAGIC 
# MAGIC - nodes_df: dataframe with three columns ('id' is the node id, 'identifier' is actual entity identifier (i.e., if it's a gene, then a gene id), and 'type' which is the entity type (e.g., 'gene'))
# MAGIC - edges_df: dataframe with three columns 'src' and 'dst' columns using unique ids from node_df, and 'type' which is the edge type (e.g., 'disease_upregulates_gene')

# COMMAND ----------

# check whether I need all these packages
from pyspark.sql import *
from pyspark.sql.functions import *
from pyspark.sql.types import *
from functools import reduce
import pandas as pd
import itertools
import copy
from pyspark.sql.functions import array, udf
from pyspark.sql.types import ArrayType, StringType, IntegerType

# COMMAND ----------

# flag for using tiny dataset
use_test_set = False

storage_account_name = "zenegraph"
storage_account_access_key = dbutils.secrets.get(scope = "zenegraph", key = "zenegraph-storage")

spark.conf.set(
  "fs.azure.account.key."+storage_account_name+".blob.core.windows.net",
  storage_account_access_key)
# needed to use pandas
spark.conf.set("spark.sql.execution.arrow.enabled", "true")

# reading the node json from blob storage
nodes_location = "wasbs://bikg-sources@zenegraph.blob.core.windows.net/hetionet_nodes.json"
nodes_df = spark.read.json(nodes_location)

# reading the edge json from blob storage
edges_location = "wasbs://bikg-sources@zenegraph.blob.core.windows.net/hetionet_edges.json"
edges_df = spark.read.json(edges_location)

# use small sample first
if use_test_set:
  nodes_df = nodes_df.limit(1000)
  edges_df = edges_df.limit(1000)

# number of nodes 
nodes_df.count()

# COMMAND ----------

# number of edges
edges_df.count()

# COMMAND ----------

# MAGIC %md
# MAGIC ## Cleaning node dataframe
# MAGIC 
# MAGIC The following are both heuristics, but will lead to having a single value per column
# MAGIC 1. From the column 'identifiers' select the shortes element of the list i.e., from ["GO:0031753","http://purl.obolibrary.org/obo/GO_0031753"] to "GO:0031753"
# MAGIC 2. From the column 'labels' select the longest element in the list i.e., from ["SERPINF2","serpin peptidase inhibitor, clade F (alpha-2 antiplasmin, pigment epithelium derived factor), member 2"] to "serpin peptidase inhibitor, clade F (alpha-2 antiplasmin, pigment epithelium derived factor)"

# COMMAND ----------

# this function takes a list of strings and returns the shortest one
def select_shortest_element(list_of_strings):
  shortest = __builtins__.min(list_of_strings, key=len)
  return shortest
# convert to udf
select_shortest_element_udf = udf(select_shortest_element, StringType())

# this function takes a list of strings and returns the longest one
def select_longest_element(list_of_strings):
  longest = __builtins__.max(list_of_strings, key=len)
  return longest
# convert to udf
select_longest_element_udf = udf(select_longest_element, StringType())

# apply functions to nodes_df for cleaning
nodes_df_cleaned = nodes_df.withColumn('identifiers', select_shortest_element_udf(nodes_df.identifiers)).withColumn('labels', select_longest_element_udf(nodes_df.labels))
display(nodes_df_cleaned)


# COMMAND ----------

# MAGIC %md
# MAGIC ## Cleaning edge dataframe
# MAGIC 
# MAGIC Rename columns as src and dst which is what GraphFrames needs.

# COMMAND ----------

edges_df_cleaned = edges_df.selectExpr('source_id as src', 'target_id as dst', 'type as type')
display(edges_df_cleaned)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Export to csv for construction of graph and further analysis

# COMMAND ----------

def write_to_azure_blob(df,filename):
  write_location = "wasbs://bikg-sources@zenegraph.blob.core.windows.net/"+filename+".csv"
  df.coalesce(1).write.format("com.databricks.spark.csv").option("header", "true").save(write_location)

write_to_azure_blob(nodes_df_cleaned,'hetio_nodes')
write_to_azure_blob(edges_df_cleaned,'hetio_edges')

# COMMAND ----------



# COMMAND ----------



# COMMAND ----------

