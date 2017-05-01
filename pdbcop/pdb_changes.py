import numpy as np
import pandas as pd 
from biopandas.pdb import PandasPdb

##########
# TO DO
# tolerance B, tolerance d filtering
# is left/right necessary? 
##########

def pdb_changes(Pdb1, Pdb2, record = 'ATOM'):
    """
    """
    if(Pdb1.df[record].ix[:,0:10].equals(Pdb2.df[record].ix[:,0:10])):
        changes = Pdb1.df[record].ix[:,0:10].copy()
        changes['b_factor_change'] = Pdb1.df[record]['b_factor'] - Pdb2.df[record]['b_factor']
        changes['b_factor_final'] = Pdb2.df[record]['b_factor']
        changes['dx'] = Pdb1.df[record]['x_coord'] - Pdb2.df[record]['x_coord']
        changes['dy'] = Pdb1.df[record]['y_coord'] - Pdb2.df[record]['y_coord']
        changes['dz'] = Pdb1.df[record]['z_coord'] - Pdb2.df[record]['z_coord']
        changes['coord_change'] = np.sqrt(changes['dx']**2+changes['dy']**2+\
                                 changes['dz']**2)
        changes['occupancy'] = Pdb2.df[record]['occupancy']
        right = None
        left = None
    else:  # TO DO test if it works with update waters option in refmac/phenix
        merged = Pdb1.df[record].ix[:,0:10].merge(Pdb2.df[record].ix[:,0:10], indicator=True, how='outer')
        right = merged[merged['_merge'] == 'right_only']
        left = merged[merged['_merge'] == 'left_only']
        atom_nu1 = list(left['atom_number'])
        atom_nu2 = list(right['atom_number'])
        tmp_pdb1 = Pdb1.df[record][~Pdb1.df[record].atom_number.isin(atom_nu1)]
        tmp_pdb2 = Pdb2.df[record][~Pdb2.df[record].atom_number.isin(atom_nu2)]
        changes = tmp_pdb1.ix[:,0:10].copy()
        changes['b_factor_change'] = Pdb1.df[record]['b_factor'] - Pdb2.df[record]['b_factor']
        changes['b_factor_final'] = Pdb2.df[record]['b_factor']
        changes['dx'] = tmp_pdb1['x_coord'] - tmp_pdb2['x_coord']
        changes['dy'] = tmp_pdb1['y_coord'] - tmp_pdb2['y_coord']
        changes['dz'] = tmp_pdb1['z_coord'] - tmp_pdb2['z_coord']
        changes['coord_change'] = np.sqrt(changes['dx']**2+changes['dy']**2+\
                                 changes['dz']**2)
        changes['occupancy'] = tmp_pdb2['occupancy']

    return changes, left, right

if __name__ == '__main__': #just for testing... 
    
    pf1 = PandasPdb().read_pdb('../test/data/1.pdb')
    pf2 = PandasPdb().read_pdb('../test/data/2.pdb')
    changes, left, right = pdb_changes(pf1,pf2,record="HETATM")
    filter = changes['coord_change'] >= 0.2 #0.2
    print(changes[filter].sort_values('coord_change', ascending = False))
    filter = changes['b_factor_change'] >= 2 #3.5
    print(changes[filter].sort_values('coord_change', ascending = False))
