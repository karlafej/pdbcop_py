import numpy as np
from biopandas.pdb import PandasPdb

###########
# Residues
###########

def residue_count(Pdb, record = 'ATOM'):
    """
    """
    tmp = Pdb.df[record].drop_duplicates(subset=['chain_id','residue_number'])
    chains = tmp.chain_id.unique()
    res_count = {}
    for chain_id in chains:
        res_count[chain_id] = len(tmp[tmp.chain_id == chain_id].index)
    return res_count 

def nonstandard(Pdb, record = 'ATOM'):
    """
    """
    residues_std = ('ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
                'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',
                'HOH')
    return Pdb.df[record][~Pdb.df[record].residue_name.isin(residues_std)]

###########
# Alternate conformations
###########



#############
# B - factors
#############

def bvalues(Pdb, record = 'ATOM'):
    """
    """
    B_min = Pdb.df[record]['b_factor'].min()
    B_max = Pdb.df[record]['b_factor'].max()
    B_ave = Pdb.df[record]['b_factor'].mean()
    return B_min, B_max, B_ave

def blargest(Pdb, n, record = 'ATOM'):
    """
    """
    return Pdb.df[record].nlargest(n, 'b_factor')

def blowest(Pdb, n, record = 'ATOM'):
    """
    """
    return Pdb.df[record].nsmallest(n, 'b_factor')

#########
# Occupancy
#########

def occupancy_partial(Pdb):
    occ_part = {} 
    for record in ('ATOM', 'HETATM'): #TO DO test if record exists
        filter = (Pdb.df[record]['occupancy'] != 1) & (Pdb.df[record]['occupancy'] != 0)
        occ_part[record] = Pdb.df[record][filter]
    return occ_part

def occupancy_low(Pdb, lim = 0.1):
    occ_low = {}
    for record in ('ATOM', 'HETATM'): #TO DO test if record exists
        filter = ppdb.df[record]['occupancy'] <= float(lim)
        occ_low[record] = ppdb.df[record][filter]
    return occ_low


##########
# Sequence
##########

def fasta(Pdb): 
    """
    TO DO: write to a file
    """
    sequence = Pdb.amino3to1()
    sequence_list = list(sequence.loc[sequence['chain_id'] == 'A', 'residue_name'])
    for chain_id in sequence['chain_id'].unique():
        print('\n> Chain ID: %s' % chain_id)
        print(''.join(sequence.loc[sequence['chain_id'] == chain_id, 'residue_name']))
    return None

##########
# Anisotropic ADPs
##########

def eigenvalues(Pdb, record = 'ANISOU'):
    """Thermal ellipsoid eigenvalues

    Cardano's analytical algorithm:
    http://www.apc.univ-paris7.fr/Downloads/antares/Joao/OscProb_v0/doxfiles/html/zheevc3_8cxx_source.html

    input: biopandas pdb object
    output: pandas dataframe with e1, e2, e3 values

    usage: ppdb = PandasPdb().fetch_pdb('3nir')
           e123 = eigenvalues(ppdb)

    ######
    seems faster than:
    def eigenvalues(U11,U22,U33,U12,U13,U23):
        a = np.array([[U11,U12,U13],[U12,U22,U23],[U13,U23,U33]])
        evalues = list(np.linalg.eigvalsh(a)/10000)
        return evalues  

    def eigen_df(AnisouDf, record ='ANISOU'): ### slow
        EigenDf = AnisouDf.df[record].copy()
        EigenDf[['e1', 'e2','e3']] = EigenDf.apply(lambda row: eigenvalues(row['U(1,1)'], row['U(2,2)'], \
                                 row['U(3,3)'], row['U(1,2)'], \
                                 row['U(1,3)'], row['U(2,3)']), axis = 1).apply(pd.Series)
        EigenDf.drop(['record_name', 'blank_1', 'blank_2', 'blank_3', 'blank_4', 'charge'], inplace=True,axis=1)
        return EigenDf

    """
    sqrt_3 = np.sqrt(3)
    EigenDf = Pdb.df[record].copy()
    EigenDf['de'] = EigenDf['U(1,2)'] * EigenDf['U(2,3)'] 
    EigenDf['dd'] = EigenDf['U(1,2)']**2 
    EigenDf['ee'] = EigenDf['U(2,3)']**2 
    EigenDf['ff'] = EigenDf['U(1,3)']**2 
    EigenDf['m'] = EigenDf['U(1,1)'] + EigenDf['U(2,2)'] + EigenDf['U(3,3)'] 
    EigenDf['c1'] = EigenDf['U(1,1)'] * EigenDf['U(2,2)'] + \
                    EigenDf['U(1,1)'] * EigenDf['U(3,3)'] + \
                    EigenDf['U(2,2)'] * EigenDf['U(3,3)']  - \
                   (EigenDf['dd'] + EigenDf['ee'] + EigenDf['ff'] )
    EigenDf['c0'] = EigenDf['U(3,3)'] * EigenDf['dd'] + \
                    EigenDf['U(1,1)'] * EigenDf['ee'] + \
                    EigenDf['U(2,2)'] * EigenDf['ff'] - \
                    EigenDf['U(1,1)'] * EigenDf['U(2,2)'] * EigenDf['U(3,3)'] - \
                    2.0 * EigenDf['U(1,3)'] * EigenDf['de']
    EigenDf['p'] = (EigenDf['m'] ** 2) - 3.0 * EigenDf['c1'] 
    EigenDf['q'] = EigenDf['m'] * (EigenDf['p'] - 3/2 * EigenDf['c1'] ) - \
                   27/2 * EigenDf['c0'] 
    EigenDf['sqrt_p']  = np.sqrt(abs(EigenDf['p']))
    EigenDf['phi'] = 27 * (0.25 * EigenDf['c1']**2 * \
                   (EigenDf['p'] - EigenDf['c1']) + \
                   EigenDf['c0'] * \
                  (EigenDf['q'] + 27/4 * EigenDf['c0'] ))
    EigenDf['phi'] = 1/3 * np.arctan2(np.sqrt(abs(EigenDf['phi'] )), EigenDf['q'] )
    EigenDf['c'] = EigenDf['sqrt_p'] * np.cos(EigenDf['phi'] )
    EigenDf['s'] = 1/sqrt_3 * EigenDf['sqrt_p'] * np.sin(EigenDf['phi'] )

    EigenDf['e2'] = 1/3 * (EigenDf['m'] - EigenDf['c'])
    EigenDf['e1'] = (EigenDf['e2'] + EigenDf['c'])/10000
    EigenDf['e3'] = (EigenDf['e2'] + EigenDf['s'])/10000
    EigenDf['e2'] = (EigenDf['e2'] - EigenDf['s'])/10000
    EigenDf.drop(['de', 'dd', 'ee', 'ff', 'm', 'c1', 'c0', 'p', 'q', 'phi', 'c', 's','sqrt_p'],inplace=True,axis=1)
    EigenDf.drop(['record_name', 'blank_1', 'blank_2', 'blank_3', 'blank_4', 'charge'], inplace=True,axis=1)
    return EigenDf


if __name__ == '__main__': #just for testing... 

    ppdb = PandasPdb().fetch_pdb('3mtn')

    print(blargest(ppdb, 5, record = 'HETATM'))
    print(blowest(ppdb, 5))
    print(bvalues(ppdb))
    fasta(ppdb)
    print(residue_count(ppdb))
    print(nonstandard(ppdb, record = 'HETATM').head())
    print(occupancy_partial(ppdb)['ATOM'].head())
    print(occupancy_low(ppdb, 0.16))

    ppdb2 = PandasPdb().fetch_pdb('3nir')
    e123 = eigenvalues(ppdb2)
    print(e123.head())