# This is a script to download the covariance and correlation matrices published
# by the FYCoM group https://nucleardata.berkeley.edu/FYCoM/ owned by Eric F. Matthews
# These matrices are prepackaged in the code in conflux/fissionDB

import sys
import urllib.request
import numpy
import os.path

urlBase = 'https://raw.githubusercontent.com/efmatthews/FYCoM/master/matrices/'
dataCat = ['ENDF/', 'JEFF/']
dataType = 'cumulative/'
covcorr = ['cov', 'corr', 'normed_cov']
energy = ['T','F','H','SF']
Z = {
    96:'Cm',
    97:'Bk',
    98:'Cf',
    99:'Es',
    100:'Fm',
}

if __name__ == "__main__":
    targetDir = sys.argv[1]
    for i in range(220, 260):
        for z in Z:
            for e in energy:
                for mat in covcorr:
                    for cat in dataCat:
                        urlName = urlBase+cat+dataType+Z[z]+str(i)+e+'_cml_'+mat+'.csv'
                        localName=targetDir+cat+mat+"_nfy_"+str(z)+"_"+Z[z]+"_"+str(i)+"_"+e+'.csv'
                        print("Downloading: "+urlName)
                        print("Saving to: "+localName)
                        if os.path.isfile(localName):
                            continue
                        try:
                            crawl_response = urllib.request.urlopen(urlName, timeout = 30)
                            
                            urllib.request.urlretrieve(urlName, localName)
                        except urllib.error.HTTPError as err:
                            if err.code == 404:
                                print('File does not exist')
                            
