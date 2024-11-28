from wcvpy.wcvp_download import get_all_taxa, wcvp_accepted_columns, plot_native_number_accepted_taxa_in_regions

if __name__ == '__main__':
    new_taxa = get_all_taxa(accepted=True)
    to_plot = new_taxa[~new_taxa[wcvp_accepted_columns['species']].isna()]
    plot_native_number_accepted_taxa_in_regions(to_plot, wcvp_accepted_columns['species'],
                                                '', 'all_species_native_distribution.jpg',
                                                include_extinct=True)
