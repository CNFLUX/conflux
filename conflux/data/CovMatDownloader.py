# This is a script to download the covariance and correlation matrices published
# by the FYCoM group https://nucleardata.berkeley.edu/FYCoM/ owned by Eric F. Matthews
# These matrices are prepackaged in the code in conflux/fissionDB

import sys
import urllib.request
import numpy
import os.path
import os

urlBase = 'https://raw.githubusercontent.com/efmatthews/FYCoM/master/matrices/'
DB_path = os.environ["CONFLUX_DB"]+"/fissionDB/"
dataCat = ['/ENDF/', '/JEFF/']
dataType = 'cumulative/'
covcorr = ['cov', 'corr', 'normed_cov']
energy = ['T','F','H','SF']
Z = {
    90:'Th',
    92:'U',
    94:'Pu'
    # 96:'Cm',
    # 97:'Bk',
    # 98:'Cf',
    # 99:'Es',
    # 100:'Fm',
}

if __name__ == "__main__":
    for z in Z:
        for i in range(z*2+47, z*2+56):
            for e in energy:
                for mat in covcorr:
                    for cat in dataCat:
                        urlName = urlBase+cat+dataType+Z[z]+str(i)+e+'_cml_'+mat+'.csv'
                        localName=DB_path+cat+mat+"_nfy_"+str(z)+"_"+Z[z]+"_"+str(i)+"_"+e+'.csv'
                        if os.path.isfile(localName):
                            print('File exists: '+localName)
                            continue
                        try:
                            crawl_response = urllib.request.urlopen(urlName, timeout = 30)
                            print("Downloading: "+urlName)
                            urllib.request.urlretrieve(urlName, localName)
                            print("Saving to: "+localName)
                        except urllib.error.HTTPError as err:
                            if err.code == 404:
                                print(cat+dataType+Z[z]+str(i)+e+'_cml_'+mat+' Matrices do not exist')
