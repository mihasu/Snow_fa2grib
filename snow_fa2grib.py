#!/lustre/home/hirlam2/.conda/envs/epygram/bin/python3

import os
import subprocess
import epygram
import numpy as np
import warnings
from scipy.interpolate import griddata
import sys

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

# Define file paths and directories
input_directory = '/lustre/tmp/hirlam2/Snowparams_BMM/metcoop_input/' 
climate_file_path = '/lustre/tmp/cooper/harmonie/MEPS_prod/climate/METCOOP25D/Const.Clim.sfx'
reference_grib2_path = '/lustre/tmp/hirlam2/Snowparams_BMM/metcoop_input/reference_grib2'

# Define grib keys for new parameters
fid_snow_depth = {
    'FA': 'SNOW_DEPTH',
    'generic': {
        'discipline': 0,
        'parameterCategory': 1,
        'parameterNumber': 11,
        'typeOfFirstFixedSurface': 1,
    }
}
fid_snow_fraction = {
    'FA': 'SNOW_FRACTION',
    'generic': {
        'discipline': 0,
        'parameterCategory': 1,
        'parameterNumber': 42,
    }
}
fid_ice_fraction = {
    'FA': 'ICE_FRACTION',
    'generic': {
        'discipline': 10,
        'parameterCategory': 2,
        'parameterNumber': 0,
        'typeOfFirstFixedSurface': 1,
    }
}

def calculate_fields(fields):
    snow_depth = np.zeros_like(fields['h_snow'])
    snow_fraction = np.zeros_like(fields['h_snow'])
    ice_fraction = np.zeros_like(fields['h_ice'])
    inland_ice_frac = np.zeros_like(fields['h_ice'])

    cond1 = (fields['frac_town'] + fields['frac_sea']) == 1.0
    cond2 = (fields['frac_town'] < 1.0) & (fields['frac_nature'] == 0.0) & (fields['frac_water'] > 0.0)
    cond3 = (fields['frac_town'] < 1.0) & (fields['frac_nature'] > 0.0)

    snow_depth[cond1] = 0.0
    snow_depth[cond2] = fields['h_snow'][cond2]
    snow_depth[cond3] = (
        (fields['h_snow'][cond3] * fields['frac_water'][cond3]) +
        (fields['dsn_t_isba'][cond3] * (fields['frac_nature'][cond3] + fields['frac_town'][cond3]))
    ) / (fields['frac_water'][cond3] + fields['frac_nature'][cond3] + fields['frac_town'][cond3])

    snow_fraction[cond1] = 0.0
    cond4 = cond2 & (fields['h_snow'] == 0.0)
    cond5 = cond2 & (fields['h_snow'] > 0.0)
    cond6 = cond3 & (fields['h_snow'] == 0.0)
    cond7 = cond3 & (fields['h_snow'] > 0.0)

    snow_fraction[cond4] = 0.0
    snow_fraction[cond5] = 1.0
    snow_fraction[cond6] = fields['psn_isba'][cond6] * (fields['frac_nature'][cond6] + fields['frac_town'][cond6])
    snow_fraction[cond7] = (
        fields['psn_isba'][cond7] * (fields['frac_nature'][cond7] + fields['frac_town'][cond7]) +
        fields['frac_water'][cond7]
    )

    cond8 = fields['frac_water'] > 0.0
    cond9 = cond8 & (fields['h_ice'] == 0.0)
    cond10 = ((fields['frac_sea'] + fields['frac_water']) > 0.0) & ((fields['frac_sea'] + fields['frac_water']) <= 1.0)
    
    inland_ice_frac[cond9] = 1.0
    inland_ice_frac[cond8 & ~cond9] = 0.0
    ice_fraction[cond10] = (
        fields['sic'][cond10] * fields['frac_sea'][cond10] +
        inland_ice_frac[cond10] * fields['frac_water'][cond10]
    )
    return snow_depth, snow_fraction, ice_fraction

def write_field(output_fa, fa, field_name, data, fid_info):
    field = fa.readfield(field_name)
    field.setdata(data)
    field.fid = fid_info
    output_fa.writefield(field)

def resample_field(field_data, source_geometry, target_geometry):
    src_lon, src_lat = source_geometry.get_lonlat_grid()
    tgt_lon, tgt_lat = target_geometry.get_lonlat_grid()

    points = np.column_stack((src_lon.ravel(), src_lat.ravel()))
    values = field_data.ravel()
    target_points = np.column_stack((tgt_lon.ravel(), tgt_lat.ravel()))
    resampled_data = griddata(points, values, target_points, method='linear')
    resampled_data = resampled_data.reshape(tgt_lon.shape)
    return np.nan_to_num(resampled_data, nan=0.0)

def write_output_files(output_fa_path, output_grib2_path, output_grib2_tmp, fields, snow_depth, snow_fraction, ice_fraction, ref_geometry, fa):
    with epygram.formats.resource(output_fa_path, 'w', fmt='FA') as output_fa:
        output_fa.geometry = fields['geometry']
        output_fa.validity = fields['validity']

        write_field(output_fa, fa, 'SFX.H_SNOW', snow_depth, {'FA': 'SNOW_DEPTH'})
        write_field(output_fa, fa, 'SFX.PSN_ISBA', snow_fraction, {'FA': 'SNOW_FRACTION'})
        write_field(output_fa, fa, 'SFX.SIC', ice_fraction, {'FA': 'ICE_FRACTION'})

    with epygram.formats.resource(output_fa_path, 'r', fmt='FA') as fa_file:
        with epygram.formats.resource(output_grib2_tmp, 'w', fmt='GRIB') as grib2_file:
            for fieldname in fa_file.listfields():
                field = fa_file.readfield(fieldname)
                resampled_data = resample_field(field.getdata(), field.geometry, ref_geometry)

                if fieldname == 'SNOW_DEPTH':
                    new_fid = fid_snow_depth
                elif fieldname == 'SNOW_FRACTION':
                    new_fid = fid_snow_fraction
                elif fieldname == 'ICE_FRACTION':
                    new_fid = fid_ice_fraction 

                new_field = epygram.fields.H2DField(
                    fid=new_fid,
                    geometry=ref_geometry,
                    validity=field.validity,
                    structure='H2D'
                )
                
                new_field.setdata(resampled_data)
                grib2_file.writefield(new_field, packing={'packingType': 'grid_ccsds'}, ordering={'jScansPositively': 1})
    result = subprocess.run(['grib_set', '-s',
                             'latitudeOfFirstGridPointInDegrees=50.319616,'
                             'longitudeOfFirstGridPointInDegrees=0.278280,'
                             'Latin2=63300000,'
                             'latitudeOfSouthernPoleInDegrees=-90,'
                             'NV=132,'
                             'productDefinitionTemplateNumber=1,'
                             'centre=251,'
                             'subCentre=255,'
                             'generatingProcessIdentifier=0,',
                             output_grib2_tmp, 
                             output_grib2_path],
                            capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Error running grib_set: {result.stderr}")
    else:
        print(f"grib_set command executed successfully on {output_grib2_path}")
    
    os.remove(output_fa_path)
    os.remove(output_grib2_tmp)
    print(f"Removed temporary FA file: {output_fa_path}")
    print(f"Removed temporary GRIB2 file: {output_grib2_tmp}")

def read_climate_fields(climate_file_path):
    print("Reading climate data")
    with epygram.formats.resource(climate_file_path, 'r') as climate:
        climate_fields = {
            'frac_town': climate.readfield('SFX.FRAC_TOWN').getdata(),
            'frac_sea': climate.readfield('SFX.FRAC_SEA').getdata(),
            'frac_nature': climate.readfield('SFX.FRAC_NATURE').getdata(),
            'frac_water': climate.readfield('SFX.FRAC_WATER').getdata(),
        }
    return climate_fields

def read_reference_geometry(reference_grib2_path):
    print("Reading reference GRIB2 geometry")
    with epygram.formats.resource(reference_grib2_path, 'r', fmt='GRIB') as ref_grib2:
        fields_list = ref_grib2.listfields()
        if not fields_list:
            raise ValueError("No fields found in the reference GRIB2 file.")
        reference_field = ref_grib2.readfield(fields_list[0])
        return reference_field.geometry

def process_file(fa_file_path, climate_fields, reference_geometry):
    print(f"Processing file: {fa_file_path}")

    with epygram.formats.resource(fa_file_path, 'r') as fa:
        fields = {
            'geometry': fa.geometry,
            'validity': fa.validity,
            'h_snow': fa.readfield('SFX.H_SNOW').getdata(),
            'dsn_t_isba': fa.readfield('SFX.DSN_T_ISBA').getdata(),
            'psn_isba': fa.readfield('SFX.PSN_ISBA').getdata(),
            'sic': fa.readfield('SFX.SIC').getdata(),
            'h_ice': fa.readfield('SFX.H_ICE').getdata(),
            **climate_fields
        }

        timestamp = fa.validity.get().strftime('%Y%m%d%H')
        output_fa_file_path = f'/lustre/tmp/hirlam2/Snowparams_BMM/out/snowparams_{timestamp}.fa'
        output_grib2_tmp_path = f'/lustre/tmp/hirlam2/Snowparams_BMM/out/snowparams_{timestamp}_tmp.grib2'
        output_grib2_file_path = f'/lustre/tmp/hirlam2/Snowparams_BMM/out/snowparams_{timestamp}.grib2'


    snow_depth, snow_fraction, ice_fraction = calculate_fields(fields)
    write_output_files(output_fa_file_path, output_grib2_file_path, output_grib2_tmp_path, fields, snow_depth, snow_fraction, ice_fraction, reference_geometry, fa)
    print(f"Output saved to {output_grib2_file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: snow_fa2grib.py <ICMSHSELE_file_path>")
        sys.exit(1)

    fa_file_path = sys.argv[1]
    if not os.path.isfile(fa_file_path):
        print(f"File not found: {fa_file_path}")
        sys.exit(1)

    climate_fields = read_climate_fields(climate_file_path)
    reference_geometry = read_reference_geometry(reference_grib2_path)
    process_file(fa_file_path, climate_fields, reference_geometry)
    print("Processing complete.")
