import pandas
import os


# This query represents dataset "m_dataset" for domain "person" and was generated for All of Us Controlled Tier Dataset v7
dataset_36813937_person_sql = """
    SELECT
        person.person_id,
        p_gender_concept.concept_name as gender,
        person.birth_datetime as date_of_birth,
        p_sex_at_birth_concept.concept_name as sex_at_birth 
    FROM
        `""" + os.environ["WORKSPACE_CDR"] + """.person` person 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_gender_concept 
            ON person.gender_concept_id = p_gender_concept.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_sex_at_birth_concept 
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id"""

dataset_36813937_person_df = pandas.read_gbq(
    dataset_36813937_person_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")

dataset_36813937_person_df.head(5)

print("Got dataset.")

import pandas
import os

# This query represents dataset "m_dataset" for domain "survey" and was generated for All of Us Controlled Tier Dataset v7
dataset_36813937_survey_sql = """
    SELECT
        answer.person_id,
        answer.survey_datetime,
        answer.question,
        answer.answer  
    FROM
        `""" + os.environ["WORKSPACE_CDR"] + """.ds_survey` answer   
    WHERE
        (
            question_concept_id IN (
                1585892
            )
        )"""

dataset_36813937_survey_df = pandas.read_gbq(
    dataset_36813937_survey_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")

print("Got dataset. Here is survey:")

dataset_36813937_survey_df.head(5)

m_df = dataset_36813937_survey_df.copy()

m_df.head(5)

list(m_df)

demog_df = dataset_36813937_person_df.copy()



list(demog_df)


demog_df.shape

m_demog = m_df.merge(demog_df[['person_id', 'date_of_birth', 'sex_at_birth']], on='person_id', how='left')

m_demog.head(5)

list(m_demog)

m_demog.shape

m_demog.dtypes

print(m_demog.columns)

m_demog["age"] = m_demog["survey_datetime"] - m_demog["date_of_birth"] 

m_demog.head(5)

m_demog.dtypes

import pandas as pd
m_demog.age = m_demog.age/pd.Timedelta('365 days')

m_demog.head(5)

m_demog.dtypes

m_demog.age.describe()

m_demog.hist(column = "age")

m_demog.shape

m_demog = m_demog.dropna()
m_demog.shape

m_demog.sex_at_birth.value_counts()

m_demog.loc[(m_demog['sex_at_birth'] == 'Male'), 'is_male'] = 1
m_demog.loc[~(m_demog['sex_at_birth'] == 'Male'), 'is_male'] = 0

m_demog.is_male.value_counts()

m_demog.head(5)

list(m_demog)

m_demog.shape

m_demog = m_demog[["person_id", "answer", "is_male", "age"]]

m_demog.head(5)

list(m_demog)

m_demog.shape

m_demog.to_csv("m_improved.tsv", sep = "\t")

!ls

import os
bucket = os.getenv("WORKSPACE_BUCKET")
bucket

!gsutil -m cp "m_improved.tsv" {bucket}/data/

print("Bucket:")

!gsutil -m ls {bucket}/data/

from datetime import datetime
start = datetime.now()

import os
bucket = os.getenv("WORKSPACE_BUCKET")
bucket

print("Importing hail")

import hail as hl

hl.stop()

hl.init(default_reference = "GRCh38")

mt_path = os.getenv("WGS_CLINVAR_SPLIT_HAIL_PATH")
#mt_path = os.getenv("MICROARRAY_HAIL_STORAGE_PATH")
mt_path

mt = hl.read_matrix_table(mt_path)

print("About to get short read")

#!gsutil -u $GOOGLE_PROJECT -m ls gs://fc-aou-datasets-controlled/v7/microarray/hail.mt
!gsutil -u $GOOGLE_PROJECT -m ls gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/clinvar/splitMT/hail.mt

mt.describe(widget=True)

# get number of variants and samples in the MatrixTable
mt.count()
# uncomment the command below to show first 5 rows if you have a RW account
mt.row.show(5)
# uncomment the command below to show first 5 columns if you have a RW account
mt.col.show(5)
# uncomment the command below to show GT status of the first 4 participants for first 5 variants if you have a RW account
mt.GT.show(n_rows=5, n_cols=4)

# get path of file containing list of related samples
#!gsutil -u $GOOGLE_PROJECT -m ls gs://fc-aou-datasets-controlled/v7/

# get path of file containing list of related samples
#!gsutil -u $GOOGLE_PROJECT -m ls gs://fc-aou-datasets-controlled/v7/microarray

#!gsutil -u $GOOGLE_PROJECT -m ls -r gs://fc-aou-datasets-controlled/v7/microarray/** | grep relatedness


print("Restricting to chr21...")

# restrict to small region
test_intervals = [f'chr{i}' for i in range(1, 23)]

mt = hl.filter_intervals(
    mt,
    [hl.parse_locus_interval(x,)
     for x in test_intervals])
mt.count()
# Sample QC
# get path of file containing list of related samples
!gsutil -u $GOOGLE_PROJECT -m ls gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/
# read the relatedness file using Hail
samples_to_remove_list = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv"
samples_to_remove = hl.import_table(samples_to_remove_list, key = "sample_id")
# remove samples from the MatrixTable
mt = mt.anti_join_cols(samples_to_remove)
mt.count()

print("Doing QC...")

# Variant QC
# recompute AF
mt = mt.annotate_rows(info = hl.agg.call_stats(mt.GT, mt.alleles))
# keep variants with minor allele frequency (maf) greater than 0.01
mt = mt.filter_rows(hl.min(mt.info.AF) > 0.01, keep = True)
mt.count()
# Population stratification
# get path of the gentic predicted ancestry file
!gsutil -u $GOOGLE_PROJECT -m ls gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/
# read ancestry file using Hail
ancestry_pred_path = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
ancestry_pred = hl.import_table(ancestry_pred_path,
                               key="research_id", 
                               impute=True, 
                               types={"research_id":"tstr","pca_features":hl.tarray(hl.tfloat)})

# link ancestry file to MatrixTable
mt_ancestry = mt.annotate_cols(ancestry_pred = ancestry_pred[mt.s])
# original mt
mt.col.describe()
# mt after adding ancestry info
mt_ancestry.col.describe()
# Load phenotypic data

# get path to workspace bucket
import os
bucket = os.getenv("WORKSPACE_BUCKET")
bucket
# list all files in the bucket subfolder `data`
!gsutil ls {bucket}/data/
# read pheno file
phenotype_filename = f"{bucket}/data/m_improved.tsv"

# Print every unique response in the 'answer' column
unique_responses = m_demog['answer'].unique()

print("Unique responses in the 'answer' column:")
print(unique_responses)

import matplotlib.pyplot as plt

# Count occurrences of each unique response
response_counts = m_demog['answer'].value_counts()

# Create a bar plot of the response counts
response_counts.plot(kind='bar')

# Set the title and labels
plt.title('Count of Each Unique Response in the Answer Column')
plt.xlabel('Response')
plt.ylabel('Count')

# Show the plot
plt.show()
print("Making binary...")
# Filter dataset to include only the relevant categories
filtered_demog = m_demog[
    m_demog['answer'].isin([
        'Current Marital Status: Married', 
        'Current Marital Status: Never Married', 
        'Current Marital Status: Living With Partner'
    ])
].copy()  # .copy() to explicitly create a copy

# Create the new binary variable using .loc for safe assignment
filtered_demog.loc[:, 'marital_status_binary'] = filtered_demog['answer'].apply(
    lambda x: 1 if x == 'Current Marital Status: Married' else 0
)

filtered_demog['is_male'] = filtered_demog['is_male'].astype(int)  # Convert 'is_male' to int
filtered_demog['age'] = filtered_demog['age'].astype(float)  # Ensure 'age' is float

# Drop unnecessary columns
filtered_demog = filtered_demog.drop(columns=['answer'])

# Reset index after filtering
filtered_demog.reset_index(drop=True, inplace=True)

print("New dataset:")

# Display the first few rows of the new dataset to verify
print(filtered_demog.head())

# Calculate the percentage of people who are married
percentage_married = (filtered_demog['marital_status_binary'].sum() / len(filtered_demog)) * 100

print(f'Percentage of people who are married: {percentage_married:.2f}%')

import os

# Step 1: Save the filtered DataFrame to a TSV file locally
filtered_demog.to_csv("/tmp/filtered_phenotypes.tsv", sep="\t", index=False)

# Step 2: Get the path to your workspace bucket
bucket = os.getenv("WORKSPACE_BUCKET")

# Define the destination path in your bucket
destination_path = f"{bucket}/data/filtered_phenotypes.tsv"
print("Copying data...")
# Step 3: Use gsutil to upload the file to your bucket
!gsutil cp "/tmp/filtered_phenotypes.tsv" {destination_path}

# Verify the upload was successful by listing the files in the bucket subfolder `data`
!gsutil ls {bucket}/data/

# Import the phenotype data from the uploaded TSV file, ensuring correct types
bucket = os.getenv("WORKSPACE_BUCKET")

# Define the destination path in your bucket dynamically
file_name = "filtered_phenotypes.tsv"  # Name of the file
folder_name = "data"  # Folder name in the bucket where the file is stored
phenotype_path = f"{bucket}/{folder_name}/{file_name}"

print("importing phenotypes...")
phenotypes = hl.import_table(
    phenotype_path,
    types={
        'person_id': hl.tstr, 
        'marital_status_binary': hl.tint, 
        'is_male': hl.tint,  # Convert 'is_male' to int (0 or 1)
        'age': hl.tfloat64  # Ensure 'age' is imported as float64
    },
    key='person_id'
)

# Annotate the MatrixTable with the new phenotype data
mt_pheno = mt_ancestry.annotate_cols(pheno=phenotypes[mt_ancestry.s])

# Ensure 'is_male' is treated as numerical (binary) for the regression
# Convert 'is_male' to float64 if it was imported as int
mt_pheno = mt_pheno.annotate_cols(pheno=mt_pheno.pheno.annotate(
    is_male=hl.float64(mt_pheno.pheno.is_male),
    age=hl.float64(mt_pheno.pheno.age)  # Ensure 'age' is treated as float64
))
print("Making covars...")
# Set covariates for the logistic regression, now correctly as numerical types
covariates = [
    1.0,  # Intercept
    mt_pheno.pheno.is_male,
    mt_pheno.pheno.age,
    mt_pheno.ancestry_pred.pca_features[0],
    mt_pheno.ancestry_pred.pca_features[1],
    mt_pheno.ancestry_pred.pca_features[2]
]

print("Running GWAS...")
# Conduct logistic regression with 'marital_status_binary' as the outcome
logistic_reg = hl.logistic_regression_rows(
    test='wald',
    y=mt_pheno.pheno.marital_status_binary,
    x=mt_pheno.GT.n_alt_alleles(),
    covariates=covariates
)

from bokeh.io import output_notebook, show
from bokeh.plotting import output_file, save
import os

# Initialize Bokeh output to display plots within the Jupyter notebook
output_notebook()

# Visualization: Manhattan plot
p = hl.plot.manhattan(logistic_reg.p_value)
show(p)

# Visualization: QQ plot
p2 = hl.plot.qq(logistic_reg.p_value)
show(p2)
print("Plotting...")
# Define the path dynamically for saving the Manhattan plot
manhattan_plot_filename = "manhattan_marital_status.html"
manhattan_plot_path = f"/tmp/{manhattan_plot_filename}"

# Save the Manhattan plot locally
output_file(manhattan_plot_path)
save(p)

# Get the workspace bucket path dynamically
bucket = os.getenv("WORKSPACE_BUCKET")

# Define the destination path in the bucket dynamically
results_folder = "results"  # Assuming 'results' is your desired folder in the bucket
manhattan_plot_bucket_path = f"{bucket}/{results_folder}/{manhattan_plot_filename}"

# Use gsutil to upload the Manhattan plot to your bucket
!gsutil cp {manhattan_plot_path} {manhattan_plot_bucket_path}

# Verify the upload by listing files in the results folder of the bucket
!gsutil ls {bucket}/{results_folder}/

# Construct the full path of the file in the bucket dynamically
manhattan_plot_bucket_path = f"{bucket}/{results_folder}/{manhattan_plot_filename}"

# Define the local path where you want to save the file
local_manhattan_plot_path = f"./{manhattan_plot_filename}"

print("Saving html...")
# Use the gsutil command to copy the file from the bucket to the local path
!gsutil cp {manhattan_plot_bucket_path} {local_manhattan_plot_path}

# Now, click the jupyter logo and go to the files, and find manhattan_marital_status.html.

print(“Finished!”)
