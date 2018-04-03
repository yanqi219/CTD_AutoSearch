#################################################
# Comparative Toxicogenomics Database AutoSearch
# Author: Qi Yan
# Date: 03/17/2013
#################################################
def ctdSearch(input_loc, output_loc, tempfile_loc, disease, threshold, mode):

    import urllib.request
    from bs4 import BeautifulSoup
    import requests
    import pandas as pd
    import numpy as np

    # Global variables
    # output_loc = 'C:/Users/QiYan/Desktop/'
    # tempfile_loc = 'C:/Users/QiYan/Desktop/test/'
    # disease = 'Parkinson Disease'
    # threshold = 5
    # input_loc = "C:/Users/QiYan/Desktop/genelist.txt"
    if mode == "gtc":
        gene_df = pd.read_csv(input_loc, sep="\t")
        output_df = pd.DataFrame(columns=['Disease', 'Gene Name', 'Gene ID', 'Chemical Name', 'Chemical ID', 'Inference Score'])
        j = 0
        for g in range(0, len(gene_df)):
            gene_num = gene_df.loc[g]["EntrezID"]
            gene = gene_num.astype(str)
            gene_name = gene_df.loc[g]["gene"]
            # Access the ctd database, auto search the targeted genes and get the interacted chemicals list
            ctd_url = 'http://ctdbase.org/detail.go?acc='+gene+'&view=chem&page=1&type=gene'     ###### GENE is a list, for loop needed!!
            outloc = tempfile_loc+gene+'.csv'

            page = requests.get(ctd_url)
            ##page = requests.get('http://ctdbase.org/detail.go?acc=6622&view=chem&page=1&type=gene')
            ## print(page.text[:500])
            soup = BeautifulSoup(page.text, 'html.parser')
            subsection_list = soup.find_all('div', class_='exportlinks')
            ## print(len(subsection_list))
            ## print(subsection_list)
            csv = subsection_list[0]
            suburl_csv = csv.a
            suburl_csv_plain = suburl_csv.get('href')
            suburl = 'http://ctdbase.org' + suburl_csv_plain

            # Save the chemicals list as .csv file
            urllib.request.urlretrieve(suburl, outloc)

            # Read csv file into pandas
            chem_list = pd.read_csv(outloc)
            # Search inference score of chemicals and PD, same methods as above, download the csv file
            for i in range(0, len(chem_list)):
                chem_name = chem_list.loc[i]['Chemical Name']   ###### will need a fpr loop here
                chem_id = chem_list.loc[i]['Chemical ID']


                chem_url = 'http://ctdbase.org/detail.go?acc='+chem_id+'&view=disease&page=1&type=chem'
                chem_outloc = tempfile_loc+chem_id+'.csv'

                chem_page = requests.get(chem_url)
                chem_soup = BeautifulSoup(chem_page.text, 'html.parser')
                chem_subsection_list = chem_soup.find_all('div', class_='exportlinks')
                chem_csv = chem_subsection_list[0]
                chem_suburl_csv = chem_csv.a
                chem_suburl_csv_plain = chem_suburl_csv.get('href')
                chem_suburl = 'http://ctdbase.org' + chem_suburl_csv_plain

                # Save the chemicals list as .csv file
                urllib.request.urlretrieve(chem_suburl, chem_outloc)
                # read the csv file and retrieve the inference score
                score_list = pd.read_csv(chem_outloc)
                score_pd = score_list.loc[score_list['Disease Name'] == disease]['Inference Score']
                highest_score = score_pd.values[0]
                if highest_score >= threshold:
                    output_df.loc[j] = [disease, gene_name, gene, chem_name, chem_id, highest_score]
                    j = j+1

        output_df.to_csv(output_loc+'related chemicals.csv', sep=',', index=False)
    elif mode == "ctg":
        chem_df = pd.read_csv(input_loc, sep="\t")
        output_df = pd.DataFrame(columns=['Disease', 'Chemical Name', 'Chemical ID', 'Gene Name', 'Gene ID', 'Inference Score'])
        j = 0
        for g in range(0, len(chem_df)):
            chem_num = chem_df.loc[g]["MeSHID"]
            #chem = chem_num.astype(str)
            chem_name = chem_df.loc[g]["chemical"]
            ctd_url = 'http://ctdbase.org/detail.go?type=chem&acc='+chem_num+'&view=gene'
            outloc = tempfile_loc+chem_num+'.csv'

            page = requests.get(ctd_url)
            soup = BeautifulSoup(page.text, 'html.parser')
            subsection_list = soup.find_all('div', class_='exportlinks')
            csv = subsection_list[0]
            suburl_csv = csv.a
            suburl_csv_plain = suburl_csv.get('href')
            suburl = 'http://ctdbase.org' + suburl_csv_plain

            urllib.request.urlretrieve(suburl, outloc)

            # Read csv file into pandas
            gene_list = pd.read_csv(outloc)
            # Search inference score of chemicals and PD, same methods as above, download the csv file
            for i in range(0, len(gene_list)):
                gene_name = gene_list.loc[i]['Gene Symbol']   ###### will need a fpr loop here
                gene_id = gene_list.loc[i]['Gene ID']
                gene = gene_id.astype(str)

                gene_url = 'http://ctdbase.org/detail.go?type=gene&acc='+gene+'&view=disease'
                gene_outloc = tempfile_loc+gene_name+'.csv'

                gene_page = requests.get(gene_url)
                gene_soup = BeautifulSoup(gene_page.text, 'html.parser')
                gene_subsection_list = gene_soup.find_all('div', class_='exportlinks')
                gene_csv = gene_subsection_list[0]
                gene_suburl_csv = gene_csv.a
                gene_suburl_csv_plain = gene_suburl_csv.get('href')
                gene_suburl = 'http://ctdbase.org' + gene_suburl_csv_plain

                # Save the chemicals list as .csv file
                urllib.request.urlretrieve(gene_suburl, gene_outloc)
                # read the csv file and retrieve the inference score
                score_list = pd.read_csv(gene_outloc)
                score_pd = score_list.loc[score_list['Disease Name'] == disease]['Inference Score']
                if score_pd.empty == True:
                    continue
                highest_score = score_pd.values[0]
                if type(highest_score) == type(str()):
                    highest_score = float(highest_score)
                if highest_score >= threshold:
                    output_df.loc[j] = [disease, chem_name, chem_num, gene_name, gene_id, highest_score]
                    j = j+1

        output_df.to_csv(output_loc+'related genes.csv', sep=',', index=False)
