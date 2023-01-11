"""

Code used in the calculation of single chemical hazard scores and mixtures of chemicals
for the DTT superfund project.

"""
import pandas as pd
import numpy as np
import re
from pathlib import Path
from casregnum import CAS


def count_organ_sites(sites):
    """ count how many sites are present """
    if sites==None or sites=="":
        return 0
    sstoks = sites.split(",")
    #print("Orgas {}, len {}".format(sites,len(sstoks)))
    return len(sstoks)


def gen_hazard_level(row):
    """
    :param row: Compute the hazard level for a given row in the sf_hazard_rcgia.csv data frame
    :return: level

    See: ./docs/HazardSingleChemical.xlsx for level assignment rules

    """
    # check ROC evidence
    roc = row['roc_15th_listing']
    if ("Known" in roc):
        return 1
    elif ("RAHC" in roc):
        return 2

    # get cebs genotox evidence
    if (row['cebsgt_flag'] == "Evidence"):
        cebsgt = 2
    elif (row['cebsgt_flag'] == 'Weak Evidence'):
        cebsgt = 1
    elif (row['cebsgt_flag'] == "Inadequate") or (row['cebsgt_flag'] == "No Evidence") or (row['cebsgt_flag'] == "Missing"):
        cebsgt = 0
    else:
        cebsgt = 0

    # now get organ info
    cnt_list = [0,0,0,0]
    cnt_list[0] = count_organ_sites(str(row['cebs_male_rats']))
    cnt_list[1] = count_organ_sites(str(row['cebs_female_rats']))
    cnt_list[2] = count_organ_sites(str(row['cebs_male_mice']))
    cnt_list[3] = count_organ_sites(str(row['cebs_female_mice']))
    num_sex_species = len([i for i in cnt_list if i > 0])
    multisite = any(i >= 2 for i in cnt_list)
    multisex_samespecies = (cnt_list[0]>=1 and cnt_list[1]>=1) or (cnt_list[2]>=1 and cnt_list[3]>=1)

    # now apply rules
    if cebsgt == 0:
        if not multisite:
            if num_sex_species==0:
                return 5
            elif num_sex_species==1:
                return 4
            elif num_sex_species>=2:
                return 3
        elif multisite:
            return 2
    elif cebsgt == 1:
        if not multisite:
            if num_sex_species==0 or num_sex_species==1:
                return 4
            elif num_sex_species==2 and multisex_samespecies==True:
                return 3
            elif num_sex_species==2 and multisex_samespecies==False:
                return 2
            elif num_sex_species>=3:
                return 2
        elif multisite:
            if multisex_samespecies:
                return 1
            elif not multisex_samespecies:
                return 2
    elif cebsgt == 2:
        if not multisite:
            if num_sex_species==0 or num_sex_species==1:
                return 3
            elif num_sex_species==2 and multisex_samespecies==True:
                return 3
            elif num_sex_species==2 and multisex_samespecies==False:
                return 3
            elif num_sex_species>=3:
                return 2
        elif multisite:
            if multisex_samespecies:
                return 1
            elif not multisex_samespecies:
                return 2
    print("\n\n\nNO LEVEL FOUND, row is {}".format(row))
    return -1

def compute_hazard_level():
    """ computer per chemical hazard levels
    """
    df = pd.read_csv(Path("./sf_hazard_riagc.csv"),dtype=str,na_filter=False)
    df['hazard_level'] = df.apply(gen_hazard_level, axis=1)
    df.to_csv(Path("./sf_hazard_riagch.csv"), index=False)
    print("Finished")


compute_hazard_level()