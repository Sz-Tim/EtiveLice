

import os
import xarray as xr
import argparse


#%%
def main(args):
    # Define the source and destination directories
    source_dir = f'{args.source_dir}/netcdf_{args.year}'
    destination_dir = f'{args.dest_dir}/netcdf_{args.year}'

    # Create the destination directory if it doesn't exist
    os.makedirs(destination_dir, exist_ok=True)

    # List of variables to copy
    variables_to_copy = ['u', 'v', 'ww', 'temp', 'salinity', 'short_wave', 'kh', 'viscofh']

    # Iterate over each file in the source directory
    filenames = os.listdir(source_dir)
    filenames.sort(reverse=args.rev)
    for filename in filenames:
        if filename.endswith('.nc') and not os.path.exists(os.path.join(destination_dir, filename)):
            print(f'Starting {filename}')
            # Open the source netCDF file using xarray
            src_path = os.path.join(source_dir, filename)
            ds = xr.open_dataset(src_path, decode_times={'Itime2': False})
            # Select only the specified variables
            ds_selected = ds[variables_to_copy]
            # Save the selected variables to a new netCDF file in the destination directory
            dest_path = os.path.join(destination_dir, filename)
            ds_selected.to_netcdf(dest_path)
            print(f'  Copied {filename} to {dest_path}')


#%%
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', type=int, help='WeStCOMS year', default=2024)
    parser.add_argument('--dest_dir', type=str, help='Destination directory', default='./WeStCOMS2/Archive')
    parser.add_argument('--source_dir', type=str, help='Source directory', default='/media/archiver/common/sa01da-work/WeStCOMS2/Archive')
    parser.add_argument('--rev', type=bool, help='Reverse download order?', default=False)
    args = parser.parse_args()
    main(args)
