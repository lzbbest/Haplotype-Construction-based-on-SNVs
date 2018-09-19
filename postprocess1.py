# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np



cnv1 = 'No_3_Copynumbers.csv' # real segment ploidy
cnv2 = 'No_3_segments.list' # predict segment ploidy

cnv1 = pd.read_csv(cnv1,sep='\t')
cnv2 = pd.read_csv(cnv2,sep='\t')

cnv=pd.merge(cnv1,cnv2,on='start')
cnv.to_csv('No_3_compare.csv',index=False)

