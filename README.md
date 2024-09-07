# EXPRES-Tools
Tools for observing with EXPRES @ LDT
====================

Currently just a target request generator but plans to add more. 

## Target List Generator:

To generate a list with requests you will need a CSV file that defines the Star name, RA/Dec, Vmag, expected exposure times, number of exposures, number of exposure meter counts, and optional comments. An example can be found [here](https://docs.google.com/spreadsheets/d/e/2PACX-1vRBx1Q26pa51QyiDdIfm-f0kFNy1WNvIBXu73HvWoPqu7Q9luI9av4UNeQ1Id_nfTWrEQv6VHy0KltG/pubhtml?gid=1181914448&single=true "EXPRES Request Example")

Usage::

  from expres_utils import TargetList
  import pandas as pd
  targets = pd.read_csv(expres_requests.csv,comment='#')
  
  tl = TargetList(targets=targets)
  tl.make_targetlist(target_list_name='example_requests',save=True)

