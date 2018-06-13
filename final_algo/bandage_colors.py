# coding: utf-8

from read_files import *
import os


def color(dataset_name="example"):
    print("data/{}/bandage_colors".format(dataset_name))
    if not os.path.exists("data/{}/bandage_colors".format(dataset_name)):
        os.makedirs("data/{}/bandage_colors".format(dataset_name))

    G = read_graph(dataset_name)
    df_ref, df_desman, strains = read_answers(G, dataset_name)

    for cur_s in strains:

        df_ref['color'] = "#b0b0b0"  # grey

        long = df_ref['length'] >= 0 #1500
        single = df_ref['single_copy']
        real_true = df_ref[cur_s] == 1
        desman_true = df_desman[cur_s] == 1

        df_ref.loc[~long & real_true, 'color'] = 'Brown'

        df_ref.loc[long & single & real_true & desman_true, 'color'] = 'Lime'
        df_ref.loc[long & ~single & real_true & desman_true, 'color'] = 'Green'

        df_ref.loc[long & single & real_true & ~desman_true, 'color'] = 'Teal'
        df_ref.loc[long & ~single & real_true & ~desman_true, 'color'] = 'Navy'

        df_ref.loc[long & single & ~real_true & desman_true, 'color'] = 'Yellow'
        df_ref.loc[long & ~single & ~real_true & desman_true, 'color'] = 'Orange'

        df_ref['strains_print'] = df_ref['strains'].apply(lambda x: ", ".join('{}({})'.format(k, v) for k, v in x.items()))
        df_ref['strains_print'] = df_ref['strains_print'].apply(lambda x: x.replace('(1)', ''))

        #df_ref[['strains_print', 'color']].to_csv("bandage_colors/{}.csv".format(cur_s), index_label='name')

        df_color_1 = df_ref[['strains_print', 'color']]
        df_color_1.index = df_color_1.index.astype(str) + '+'

        df_color_2 = df_ref[['strains_print', 'color']]
        df_color_2.index = df_color_2.index.astype(str) + '-'

        df_color = pd.concat([df_color_1, df_color_2])

        df_color.to_csv("data/{}/bandage_colors/{}.csv".format(dataset_name, cur_s), index_label='name')


color("example")
