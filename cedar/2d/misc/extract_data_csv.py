import os
import re
import csv
from pathlib import Path

# Root directory
root = Path('/project/def-goluskin/ollie/ollie_rb_data/redoing_test')

# Regex patterns to extract data
re_ra_pr = re.compile(r"Ra = ([\d.eE+-]+), Pr = ([\d.eE+-]+)")
re_time = re.compile(r"Simulation end time: ([\d.eE+-]+)")
re_avg_start = re.compile(r"Time averaging begins at: ([\d.eE+-]+)")
re_max_diff_sec = re.compile(r"Maximum percent difference in Nusselt over (\d) sections: ([\d.eE+-]+)%")
re_nu_inst = re.compile(r"Nu as calculated by average of instantaneous Nu: ([\d.eE+-]+)")
re_nu_cum = re.compile(r"Nu as calculated by cumulative average: ([\d.eE+-]+)")


# Walk through RA folders
for ra_dir in root.glob('ra1e*'):
    data_rows = []
    
    for pr_dir in ra_dir.glob('pr*'):
        pr_value = None
        
        for subdir in pr_dir.glob('*_Gam2'):
            info_file = subdir / 'info.txt'
            if not info_file.exists():
                continue

            with open(info_file) as f:
                content = f.read()

            # Extract data
            ra_pr_match = re_ra_pr.search(content)
            time_match = re_time.search(content)
            avg_start_match = re_avg_start.search(content)
            max_diff = re_max_diff_sec.search(content)
            nu_inst_match = re_nu_inst.search(content)
            nu_cum_match = re_nu_cum.search(content)

            row = {
                'Pr': float(ra_pr_match.group(2)) if ra_pr_match else None,
                'Sim end time': float(time_match.group(1)) if time_match else None,
                'Time avg start': float(avg_start_match.group(1)) if avg_start_match else None,
                'Num Sections': float(max_diff.group(1)),
                'Max %diff': float(max_diff.group(2)),
                'Nu inst avg': float(nu_inst_match.group(1)) if nu_inst_match else None,
                'Nu cum avg': float(nu_cum_match.group(1)) if nu_cum_match else None,
            }

            data_rows.append(row)

    # Write each RA's data to a CSV
    if data_rows:
        csv_filename = ra_dir.name + '.csv'
        with open(csv_filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=data_rows[0].keys())
            writer.writeheader()
            writer.writerows(data_rows)

        print(f"Written: {csv_filename}")
