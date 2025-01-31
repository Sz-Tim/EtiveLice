#!/bin/bash

python ./copy_trimmed_FVCOM_nc.py --year 2023 --source_dir 'D:/hydroOut/WeStCOMS2/Archive' --dest_dir 'D:/hydroOut/WeStCOMS2/Archive_temp'

python ./copy_trimmed_FVCOM_nc.py --year 2023 --source_dir 'D:/hydroOut/etive28/Archive' --dest_dir 'D:/hydroOut/etive28/Archive_temp'
