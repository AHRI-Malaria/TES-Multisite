
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Patch, Rectangle
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
# Step 1: Read the dataset from the CSV file
#hap_df = pd.read_csv('mad4_hatter_hap1.csv')

day3_positivity = pd.read_csv('../data/Map_Plot/Input/day3_positivity.csv')
acpr_df = pd.read_csv('../data/Map_Plot/Input/RDT.csv')

hap_df = pd.read_csv('../data/Map_Plot/Input/tes_haplotype_report5%101925.csv')
hap_df = hap_df[hap_df['Day'] == 0]
hap_df = hap_df[hap_df['SITE'] != 'Djibouti']

snp_df = pd.read_csv('../data/Map_Plot/Input/tes_snp_report5%101925.csv')

snp_df = snp_df[snp_df['SITE'] != 'Djibouti']
hrp_df = pd.read_csv('../data/Map_Plot/Input/hrp_df101925.csv')
# Filter rows where 'column_name' contains 'D0'
hrp_df = hrp_df[hrp_df['SampleID'].str.contains("D0")]


# print(snp_df )
# hrp_df = hrp_df[hrp_df['site'] != 'district']
# hrp_df = hrp_df[hrp_df['SampleID'].str.endswith('D0')]
# hap_df = hap_df[hap_df['Day'] == 0]

# Step 2: Process the data
grouped_df = hap_df.groupby(['SITE', 'lat', 'lon', 'gene', 'haplotype', 'Type']).size().reset_index(name='Prevalence')
#grouped_df = grouped_df[grouped_df['SITE'] != 'Djibuti']
grouped_df = grouped_df[grouped_df['gene'] != 'Pfdhfr_dhps']


# SNP data
snp_df = snp_df[snp_df['Day'] == 0]
print(snp_df )
snp_df = snp_df[snp_df['Gene'] == 'k13']
print(snp_df)

snp_df = snp_df[(snp_df['AA_CHANGE'] == 'A578S') | (snp_df['AA_CHANGE'] == 'R622I') | (snp_df['AA_CHANGE'] == 'A675V')|(snp_df['AA_CHANGE'] == 'A578A') | (snp_df['AA_CHANGE'] == 'R622R') | (snp_df['AA_CHANGE'] == 'A675A')]

specified_changes = ['R622I', 'A675V', 'A578S']

# Identify which SampleID have at least one of the specified changes
mask = snp_df.groupby('SampleID')['AA_CHANGE'].transform(lambda x: any(aa in specified_changes for aa in x))

# Update AA_CHANGE to 'WT' if the SampleID doesn't have any of the specified changes
snp_df.loc[~mask, 'AA_CHANGE'] = 'WT'




# Define site order and gene order
site_order = ['Maksegnit', 'Asayita', 'Asosa', 'Mizan', 'Abobo']
genes = grouped_df['gene'].unique()
num_genes = len(genes)
num_sites = len(site_order)

# Step 3: Define colors for haplotypes and sites
haplotypes = grouped_df['haplotype'].unique()

# Color palette for haplotypes
color_palette = plt.cm.get_cmap('tab20b', len(haplotypes))
#haplotype_colors = {haplotype: color_palette(i) for i, haplotype in enumerate(haplotypes)}
haplotype_colors = {
    "CVIET": '#393b79',
    "CVMNK": '#5254a3',
    "IRNI": '#6b6ecf',
    "ICNI": '#9c9ede',
    "NRNI": '#637939',
    "NCSI": '#8ca252',
    # "NCNI": '#b5cf6b',
    # "NRSI": '#cedb9c',
    "CAKAA": '#8c6d31',
    "SGKAA": '#bd9e39',
    "SGEAA": '#e7ba52',
    "SAKAA": '#e7cb94',
    "AAKAA": '#ad494a',
    "SGEGA": '#d6616b',
    "NFSND": '#e7969c',
    "NYSND": '#7b4173',
    "YFSND": '#a55194',
    # 'AAKAA': '#ce6dbd',
}
# NFSND
# NYSND
# YFSND
# FYSND

# SGEAA
# SAKAA
# SGKAA
# SAKGA-
# SGKGA
# SGEGA
# AAKAA
# SAEAA
# "SAKGA": '#843c39',
# SAKGA
# SGKGA
# SGEGA-12
rar_haplotypes = ['YFSND', "FYSND", 'SAEAA','SAKGA','AAKAA','SGKGA', 'NRSI', 'NRNI','NCNI','AAKAA']

# Assign the same color for all "Wildtype" haplotypes
wildtype_color = '#007dbf'
for idx, row in grouped_df.iterrows():
    if row['Type'] == 'Wildtype':
        haplotype_colors[row['haplotype']] = wildtype_color

# Step 4: Group rare haplotypes (<5%) under "Others"
def assign_haplotype_color(haplotype, haplotype_percentage, haplotype_type):
    if haplotype_type != 'Wildtype' and haplotype in rar_haplotypes:
        return 'Others'
    return haplotype

grouped_df['Haplotype_Display'] = grouped_df.apply(
    lambda row: assign_haplotype_color(row['haplotype'], row['Prevalence'], row['Type']), axis=1
)

# Recalculate prevalences after grouping rare haplotypes
grouped_df = grouped_df.groupby(['SITE', 'lat', 'lon', 'gene', 'Haplotype_Display', 'Type']).agg({'Prevalence': 'sum'}).reset_index()

# Add gray color for "Others"
haplotype_colors['Others'] = '#A9A9A9'

# Step 5: Create colors for geographic points (one color per site)
site_colors = {
    "Maksegnit": '#0000ff',
    "Asayita": '#8b0000',
    "Asosa": '#ff8c00',
    "Mizan": '#9400d3',
    "Abobo": '#cdcd00',
    "Djibouti": 'gray'
}
site_color_map = {site: site_colors[site] for site in site_order}


fig = plt.figure(figsize=(18, 10))
ax_map = plt.axes([0.06, 0.1, 0.71, 0.9], projection=ccrs.PlateCarree())

# Load country geometries from shapefiles
shapefile_path = shpreader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries')
reader = shpreader.Reader(shapefile_path)

# Load second and third admin level shapefiles for Ethiopia
shapefile_path_admin2 = 'gadm41_ETH_shp/gadm41_ETH_1.shp'
shapefile_path_admin3 = 'gadm41_ETH_shp/gadm41_ETH_2.shp'
shapefile_path_admin4 = 'gadm41_ETH_shp/gadm41_ETH_3.shp'
reader_admin2 = shpreader.Reader(shapefile_path_admin2)
reader_admin3 = shpreader.Reader(shapefile_path_admin3)
reader_admin4 = shpreader.Reader(shapefile_path_admin4)

# Highlight the specified level 4 names
highlighted_names = {
    "Aysaita": '#8b0000',
    "Gonder Zuria": '#0000ff',
    "Assosa": '#ff8c00',
    "Debub Bench": '#9400d3',
    "Abobo": '#cdcd00',
    "Djibouti": 'gray'
}



# Function to highlight specific level 4 areas
def highlight_level4_areas(ax, reader, names_colors):
    for record in reader.records():
        name = record.attributes['NAME_3']  # Assuming 'NAME_4' holds level 4 names
        if name in names_colors:
            color = names_colors[name]
            ax.add_geometries([record.geometry], ccrs.PlateCarree(),
                              edgecolor=color, facecolor=color, alpha=1, linewidth=1.5)

# Color Ethiopia in light orange
def add_countryET(ax, country_name, color):
    for country in reader.records():
        if country.attributes['NAME'] == country_name:
            ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                              facecolor=color, alpha=0.7)
def add_countryDJ(ax, country_name, color):
    for country in reader.records():
        if country.attributes['NAME'] == country_name:
            ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                              facecolor=color, alpha=1)

# Add administrative boundaries for level 2
def add_admin_boundaries(ax, reader, color, alpha=0.5, linewidth=0.6):
    for record in reader.records():
        ax.add_geometries([record.geometry], ccrs.PlateCarree(),
                          edgecolor=color, facecolor='none', alpha=alpha, linewidth=linewidth)

# Add Ethiopia with a light orange fill
add_countryET(ax_map, 'Ethiopia', '#F0F8FF')
#add_countryDJ(ax_map, 'Djibouti', 'gray')

# Add admin level 2 boundaries (for context)
add_admin_boundaries(ax_map, reader_admin2, 'gray', alpha=0.7)

# Highlight the specified level 4 names
highlight_level4_areas(ax_map, reader_admin4, highlighted_names)

# Add country borders, coastlines, and other features
ax_map.add_feature(cfeature.BORDERS, linestyle='-', alpha=0.2)
ax_map.add_feature(cfeature.COASTLINE, alpha=0.1, linewidth=0.8)
ax_map.add_feature(cfeature.LAND, alpha=0.2, facecolor='lightgray')
ax_map.add_feature(cfeature.OCEAN, facecolor='#bce6f9ff', alpha=0.6)
#Aysaita,Gonder zuria, Assosa, Mizan Aman, Abobo

# Set extent for the region
ax_map.set_extent([28, 72.2, -9, 20], crs=ccrs.PlateCarree())



# Step 5 plot pie chart for day 3 positivity

day3_colors = {
    'Positive': '#ff6666',  # Red
    'Negative': '#007dbf'   # bl
}

# Ensure SITE is categorical and sorted correctly
day3_positivity['SITE'] = pd.Categorical(day3_positivity['SITE'], categories=site_order, ordered=True)
day3_positivity = day3_positivity.sort_values('SITE')

# Loop through each site to plot a separate pie chart
for i, site in enumerate(site_order):

    xpos_pie = 0.35  # Adjust X position (place before other pies)
    ypos_pie = 0.73 - i * (0.5 / len(site_order))  # Adjust Y position per site

    ax_pie = plt.axes([xpos_pie, ypos_pie, 0.07, 0.05], frameon=True)

    # Get data for the current site
    site_pie_data = day3_positivity[day3_positivity['SITE'] == site]

    if not site_pie_data.empty:
        day3_positive = site_pie_data['Day3'].values[0]  # Percent already given
        day3_negative = 100 - day3_positive  # Remaining proportion

        # Plot the pie chart for the current site
        ax_pie.pie([day3_positive, day3_negative], 
                   labels=None, 
                   colors=[day3_colors['Positive'], day3_colors['Negative']], 
                   startangle=90)  

        ax_pie.set_aspect('equal')

        # Add total count (N=) below the pie chart
        ax_map.text(xpos_pie + 0.035, ypos_pie - 0.005, f'N={int(site_pie_data["N"].values[0])}', 
                    fontsize=10, ha='center', va='top', transform=fig.transFigure)

# Label for Day 3 positivity
xpos_label = 0.37
ypos_label = 0.81  
ax_map.text(xpos_label, ypos_label, "Day 3+", fontsize=12, color='black', ha='center', transform=fig.transFigure)

# Legend for Day 3 Positivity
day3_legend_elements = [Patch(facecolor=color, label=label) for label, color in day3_colors.items()]
day3_legend_elements = [Patch(facecolor=color, label=label) for label, color in day3_colors.items()]
day3_legend_elements.insert(3, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
day3_legend_elements.insert(4, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
day3_legend_elements.insert(5, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
legend_day3 = ax_map.legend(handles=day3_legend_elements, title='Day 03 Positivity', loc='lower left', bbox_to_anchor=(0.1789, -0.005), fontsize=8, ncol=1)

# step 6 plot pie chart for ACPR


acpr_colors = {
    'Failure': '#ff6666',  
    'ACPR': '#007dbf'   
}

# Ensure SITE is categorical and sorted correctly
acpr_df['SITE'] = pd.Categorical(acpr_df['SITE'], categories=site_order, ordered=True)
acpr_df = acpr_df.sort_values('SITE')

# Loop through each site to plot a separate pie chart
for i, site in enumerate(site_order):

    xpos_pie = 0.385  # Adjust X position (place before other pies)
    ypos_pie = 0.73 - i * (0.5 / len(site_order))  # Adjust Y position per site

    ax_pie = plt.axes([xpos_pie, ypos_pie, 0.07, 0.05], frameon=True)

    # Get data for the current site
    site_pie_data = acpr_df[acpr_df['SITE'] == site]

    if not site_pie_data.empty:
        acpr= site_pie_data['ACPR'].values[0]  # Percent already given
        failure = 100 - acpr # Remaining proportion

        # Plot the pie chart for the current site
        ax_pie.pie([acpr, failure], 
                   labels=None, 
                   colors=[acpr_colors['ACPR'], acpr_colors['Failure']], 
                   startangle=90)  

        ax_pie.set_aspect('equal')

        # Add total count (N=) below the pie chart
        ax_map.text(xpos_pie + 0.035, ypos_pie - 0.005, f'N={int(site_pie_data["N"].values[0])}', 
                    fontsize=10, ha='center', va='top', transform=fig.transFigure)

# Label for ACPR positivity
xpos_label = 0.42
ypos_label = 0.81  
ax_map.text(xpos_label, ypos_label, "ACPR", fontsize=12, color='black', ha='center', transform=fig.transFigure)

# Legend for Day 3 Positivity
acpr_legend_elements = [Patch(facecolor=color, label=label) for label, color in acpr_colors.items()]
acpr_legend_elements.insert(3, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
acpr_legend_elements.insert(4, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
acpr_legend_elements.insert(5, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
legend_acpr = ax_map.legend(handles=acpr_legend_elements, title='ACPR', loc='lower left', bbox_to_anchor=(0.28, -0.005), fontsize=8, ncol=1)




# Step 7: Function to plot pie charts
def plot_pie(ax, sizes, colors):
    ax.pie(sizes, colors=colors, startangle=90)
    ax.set_aspect('equal')

grouped_df['SITE'] = pd.Categorical(grouped_df['SITE'], categories=site_order, ordered=True)
grouped_df = grouped_df.sort_values('SITE')
# Step 8: Plot the pie charts
for row_idx, site in enumerate(site_order):
    for col_idx, gene in enumerate(genes):
        site_gene_data = grouped_df[(grouped_df['SITE'] == site) & (grouped_df['gene'] == gene)]
        if not site_gene_data.empty:
            haplotypes = site_gene_data['Haplotype_Display']
            prevalences = site_gene_data['Prevalence']
            colors = [haplotype_colors[hap] for hap in haplotypes]
            
            # Calculate the total count for this site and gene
            total_count = prevalences.sum()

            xpos_pie = 0.4 + col_idx * (0.168 / num_genes) + 0.1
            ypos_pie = 0.73 - row_idx * (0.5 / num_sites)

            ax_pie = plt.axes([xpos_pie, ypos_pie, 0.07, 0.05], frameon=True)
            plot_pie(ax_pie, prevalences, colors)
           

            # Add total count (N=) under the pie chart
            ax_map.text(xpos_pie + 0.035, ypos_pie - 0.005, f'N={total_count}', 
                        fontsize=10, ha='center', va='top', transform=fig.transFigure)

# Step 9: Create a single legend for "Others" and format legend for a 4x4 grid
legend_elements = []
added_others = False
print(grouped_df['SITE'] == 'Djibouti')
# Add gene groups and haplotypes to the legend
for gene in genes:
    # Add gene label at the beginning of each gene group in the legend
    legend_elements.append(Patch(facecolor='none', edgecolor='none', label=f'--- {gene} ---'))

    total_prevalence_for_gene = grouped_df[grouped_df['gene'] == gene]['Prevalence'].sum()
    haplotypes_in_gene = grouped_df[grouped_df['gene'] == gene]['Haplotype_Display'].unique()

    for haplotype in sorted(haplotypes_in_gene):
        if haplotype == 'Others':
            # Add "Others" once, not per gene
            if not added_others:
                label = "Others (<5%)"
        else:
            haplotype_type = grouped_df[grouped_df['Haplotype_Display'] == haplotype]['Type'].unique()[0]
            if haplotype_type == 'Wildtype':
                label = f"{haplotype} (WT)"
            else:
                label = haplotype

            legend_elements.append(Patch(facecolor=haplotype_colors[haplotype], label=label))

# Step 10: Insert empty spaces after first two genes and before "Others"
legend_elements.insert(3, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
legend_elements.insert(4, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
# legend_elements.insert(3, Patch(facecolor='white', edgecolor='white', label=' '))  # Empty label after first gene
legend_elements.insert(9, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
legend_elements.insert(18, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
legend_elements.insert(19, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
legend_elements.insert(20, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
legend_elements.insert(21, Patch(facecolor='gray', edgecolor='none', label='Others(<5%)'))  # Empty label after first gene
legend_elements.insert(22, Patch(facecolor='none', edgecolor='none', label=''))  # Empty label after first gene

#delete "Others" from legend
# legend_elements.pop(7)

legend_hap = ax_map.legend(handles=legend_elements,title='Haplotypes', loc='lower center', bbox_to_anchor=(0.66, -0.005), fontsize=8, ncol=5)
legend_hap.get_frame().set_facecolor('#f7fafd')
legend_hap.get_frame().set_alpha(1)

# Step 11: Plot geographic points for each site
# for site in site_order:
#     site_data = grouped_df[grouped_df['SITE'] == site]
#     lat = site_data['lat'].values[0]
#     lon = site_data['lon'].values[0]
#     ax_map.plot(lon, lat, '.', color=site_color_map[site], transform=ccrs.PlateCarree())


for row_idx, site in enumerate(site_order):
    ypos_pie = 0.73 - row_idx * (0.5 / num_sites)
    
    # Add the site name with black color
    ax_map.text(0.32, ypos_pie + 0.02, site, fontsize=12, color=site_color_map[site], 
                ha='left', va='center', transform=fig.transFigure)


for col_idx, gene in enumerate(genes):
    xpos_label = 0.42 + col_idx * (0.175 / num_genes)+0.11 
    ypos_label = 0.81  # Position above the pie charts
    ax_map.text(xpos_label, ypos_label, gene, fontsize=12, color='black', ha='center', transform=fig.transFigure)


snp_df['SITE'] = pd.Categorical(snp_df['SITE'], categories=site_order, ordered=True)
snp_df = snp_df.sort_values('SITE')
print(snp_df)
pie_data = pd.pivot_table(snp_df, 
                          values='SampleID',  # Counting distinct haplotypes
                          index='AA_CHANGE', 
                          columns='SITE', 
                          aggfunc='nunique',  # Unique haplotypes per site and AA_CHANGE
                          fill_value=0)
unique_samples_per_site = snp_df.groupby('SITE')['SampleID'].nunique().reset_index()
print(unique_samples_per_site)
# Using a color map from Matplotlib, e.g., 'tab20' for AA_CHANGE colors


aa_change_colors = {
    "A578S": '#9c9ede',
    "R622I": '#843c39',
    "A675V": '#de9ed6',
    "WT": '#007dbf'
}
snp_order = ['A578S', 'R622I', 'A675V', 'WT']
snp_colors = {snp: aa_change_colors[snp] for snp in snp_order}
#fill na values with 0

pie_data = pie_data.reindex(snp_order)
pie_data = pie_data.fillna(0)
print(pie_data)
# Loop through each site to plot a separate pie chart
for i, site in enumerate(site_order):


    xpos_pie = 0.44  
    ypos_pie = 0.73 - i * (0.5 / 5)


    ax_pie = plt.axes([xpos_pie, ypos_pie, 0.07, 0.05], frameon=True)

   
    # Get data for the current site
    site_pie_data = pie_data[site]
    
    # Calculate total sample count per site (N)
    total_samples = unique_samples_per_site[unique_samples_per_site['SITE'] == site]['SampleID'].values[0]

    
    # Plot the pie chart for the current site and AA_CHANGE distribution
    ax_pie.pie(site_pie_data, 
               labels=None,  # Do not display labels directly on the pie chart
               colors=[snp_colors[j] for j in site_pie_data.index],  # Color by AA_CHANGE
               startangle=90)  # No autopct for percentage labels
    ax_pie.set_aspect('equal')
  
                # Add total count (N=) under the pie chart
    ax_map.text(xpos_pie + 0.035, ypos_pie - 0.005, f'N={total_samples}', 
                        fontsize=10, ha='center', va='top', transform=fig.transFigure)
    
xpos_label = 0.469
ypos_label = 0.81  # Position above the pie charts
ax_map.text(xpos_label, ypos_label, "Pfk13", fontsize=12, color='black', ha='center', transform=fig.transFigure)
aa_legend_elements = [Patch(facecolor=aa_change_colors[snp], label=snp) for snp in snp_order]
aa_legend_elements.insert(4, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
legend_k13 = ax_map.legend(handles=aa_legend_elements, title='Pfk13 SNPs', loc='lower left', bbox_to_anchor=(0.346, -0.005), fontsize=8, ncol=1)
legend_k13.get_frame().set_facecolor('#f7fafd')

# Replace site abbreviations in one step
site_mapping = {'GM': 'Gondar',
                'AA': 'Asayta', 
                'BA': 'Assosa', 
                'MM': 'Mizan',
                'GA': 'Abobo', 
                'DJB': 'Djibouti'
                 
                }
# hrp_df['SITE'] = hrp_df['site'].replace(site_mapping)

# # Assigning 'hrp' categories based on HRP2 and HRP3 values
# hrp_df.loc[(hrp_df['HRP2'] > 0.05) & (hrp_df['HRP3'] > 0.05), 'hrp'] = 'WT'
# hrp_df.loc[(hrp_df['HRP2'] < 0.05) & (hrp_df['HRP3'] < 0.05), 'hrp'] = 'HRP2/3 deleted'
# hrp_df.loc[(hrp_df['HRP2'] > 0.05) & (hrp_df['HRP3'] < 0.05), 'hrp'] = 'HRP3 deleted'
# hrp_df.loc[(hrp_df['HRP2'] < 0.05) & (hrp_df['HRP3'] > 0.05), 'hrp'] = 'HRP2 deleted'

# hrp_df.to_csv('hrp_df.csv', index=False)
# Create pivot table for the pie chart data
hrp_df['SITE'] = pd.Categorical(hrp_df['SITE'], categories=site_order, ordered=True)
hrp_df = hrp_df.sort_values('SITE')
pie_hrp_data = pd.pivot_table(hrp_df, 
                              values='SampleID', 
                              index='hrp', 
                              columns='SITE', 
                              aggfunc='nunique',  
                              fill_value=0)
print(pie_hrp_data)
# Calculate total unique samples per site
unique_hrp_samples_per_site = hrp_df.groupby('SITE')['SampleID'].nunique().reset_index()

# Color mapping for HRP categories
hrp_deletion_colors = {
    "HRP2/3 deleted": '#B22222',
    "HRP2 deleted": '#FF0000',
    "HRP3 deleted": '#FFCCCB',
    "WT": '#007dbf'
}
hrp_order = ['HRP2/3 deleted', 'HRP2 deleted', 'HRP3 deleted', 'WT']

hrp_colors = {hrp: hrp_deletion_colors[hrp] for hrp in hrp_order}

print(pie_hrp_data)

# Loop to plot a pie chart for each site

for i, site in enumerate(pie_hrp_data.columns):
    
    xpos_pie = 0.7  
    ypos_pie = 0.73 - i * (0.5 / len(pie_hrp_data.columns))  # Adjust for number of sites

    ax_pie = plt.axes([xpos_pie, ypos_pie, 0.07, 0.05], frameon=True)

    # Get data for the current site
    site_pie_data = pie_hrp_data[site]

    # Calculate total sample count per site (N)
    total_samples = unique_hrp_samples_per_site[unique_hrp_samples_per_site['SITE'] == site]['SampleID'].values[0]

    # Plot the pie chart for the current site
    ax_pie.pie(site_pie_data, 
               labels=None,  
               colors=[hrp_colors[j] for j in site_pie_data.index],  # Match colors to HRP categories
               startangle=90)
    ax_pie.set_aspect('equal')
    
    # Add total sample count under the pie chart
    ax_map.text(xpos_pie + 0.035, ypos_pie - 0.005, f'N={total_samples}', 
                fontsize=10, ha='center', va='top', transform=fig.transFigure)

# Add a label above the pie charts
xpos_label = 0.74
ypos_label = 0.81  
ax_map.text(xpos_label, ypos_label, "HRP2/3", fontsize=12, color='black', ha='center', transform=fig.transFigure)

# Add a legend for HRP categories
hrp_legend_elements = [Patch(facecolor=hrp_colors[hrp], label=hrp) for hrp in hrp_order]
hrp_legend_elements.insert(5, Patch(facecolor='none', edgecolor='none', label=' '))  # Empty label after first gene
legend_hrp2 = ax_map.legend(handles=hrp_legend_elements, title='HRP2/3 Deletion', loc='lower left', 
                            bbox_to_anchor=(0.8957, -0.005), fontsize=8, ncol=1)
legend_hrp2.get_frame().set_facecolor('#f7fafd')
# Add the haplotype legend to the plot
ax_map.add_artist(legend_hap)
ax_map.add_artist(legend_k13)
ax_map.add_artist(legend_hrp2)
ax_map.add_artist(legend_acpr)
ax_map.add_artist(legend_day3)



# Add gridlines with labels
gl = ax_map.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.25, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 10, 'color': 'black'}
gl.ylabel_style = {'size': 10, 'color': 'black'}


# Save or show the plot
#plt.savefig('haplotype_k13SNP_hrp_prevalence_day_0.png', dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.savefig('../data/Map_Plot/Results/resistance_markers5%11325.pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.show()
