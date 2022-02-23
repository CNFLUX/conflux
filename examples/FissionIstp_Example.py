from conflux.BetaEngine import BetaEngine
from conflux.FPYEngine import FissionIstp, FissionModel
import csv

if __name__ == "__main__":


    #Here, I've initialized U-235. Note that my atomic number
    #Must come before my mass number
    U235 = FissionIstp(92,235)

    #After I've loaded up the Isotope that I want to look at, there's a couple of
    #Things that I can do with it. Let's try and calculate the covariance matrix for
    #U235.

    #Here, I am initializing the fission database data for U235 and storing it into
    #The cumulative Fission Yield (CFPY) and Independant Fission Yield (IFPY) dictionaries
    #That are associated with the U235 FissionIstp object. Do note as well, that with both
    #LoadFissionDB, LoadCovariance, and LoadCorrelation, you need not use the pre-packaged
    #Databases that we've provided. If you have a database of your own, simply pass the
    #location of the database as a string argument for the function.


    #U235.LoaddFissionDB(DBpath = "path/to/Database")
    U235.LoadFissionDB()

    #Once I've loaded up my fission data, I'll also go ahead and load up the Covariance and
    #Correlation data. Note, that to load up the Correlation and Cumulative Fission data, I
    #Must call LoadFissionDB() first, otherwise loading up these two data files will not work.
    #Also note that the information from both of these Databases is going to be stored inside
    #Your CFPY dictionary, with an associated 'key' for the data.


    #U235.LoadCovariance(DBpath = "path/to/Database")
    U235.LoadCovariance()
    #U235.LoadCorrelation(DBpath = "path/to/Database")
    U235.LoadCorrelation()

    #Now that I've gotten all my data, it's time to calculate the covariance matrix. This is
    #Saved in the CFPY dictionary, and I can then use that specific dictionary to
    #manipulate/look at my covariance matrix. Note that the argument that I pass into
    #CalcCovariance is the neutron energy, with 0 = thermal, 5 = fast, and 14 = ultra-relativistic
    #Neutrons


    #U235.CalcCovariance(Ei=5)
    U235.CalcCovariance(Ei=0)


    #This last little block of code writes the covariance matrix out to a csv file that can
    #Then be manipulated later.
    with open('cov_235_U_processed.csv', 'w', newline='') as csvoutput:
        fieldnames = ['']
        #Note here the fact that we're looking at CFPY[0]. If we had changed our
        #Neutron energy to say 5, then it would become CFPY[5]
        for i in U235.CFPY[0]:
            fieldnames.append(i)
        writer = csv.DictWriter(csvoutput, fieldnames=fieldnames)
        writer.writeheader()
        for i in U235.CFPY[0]:
            rowcontent = U235.CFPY[0][i].cov
            rowcontent[''] = i
            writer.writerow(rowcontent)

    print("I'm donw writing out the covariance matrix!")