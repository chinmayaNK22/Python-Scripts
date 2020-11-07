import requests ## python -m pip install requests

#### enrichment ####
import json

string_api_url = "https://string-db.org/api"
output_format = "json"
method = "enrichment"


##
## Construct the request
##

request_url = "/".join([string_api_url, output_format, method])

##
## Set parameters
##

my_genes = ['7227.FBpp0074373', '7227.FBpp0077451', '7227.FBpp0077788',
            '7227.FBpp0078993', '7227.FBpp0079060', '7227.FBpp0079448']

params = {

    "identifiers" : "%0d".join(my_genes), # your protein
    "species" : 7227, # species NCBI identifier 
    "caller_identity" : "www.awesome_app.org" # your app name

}

##
## Call STRING
##

response = requests.post(request_url, data=params)

##
## Read and parse the results
##

data = json.loads(response.text)
write1 = open('test2.txt', 'w')
for row in data:
    term = row["term"]
    preferred_names = ",".join(row["preferredNames"])
    fdr = float(row["fdr"])
    description = row["description"]
    category = row["category"]
    number_of_genes = str(row["number_of_genes"])
    number_of_genes_in_background = str(row["number_of_genes_in_background"])
    inputGenes = str(row["inputGenes"])
    
    if fdr < 0.05:

        ## print significant GO Process annotations

        write1.write("\t".join([category, term, preferred_names, str(fdr), description, number_of_genes, number_of_genes_in_background, inputGenes]) + '\n')

write1.close()
