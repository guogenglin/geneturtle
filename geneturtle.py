# -*- coding: utf-8 -*-
"""
Created on Wed May 28 15:59:29 2025

@author: Genglin Guo
@e-mail: 2019207025.njau.edu.cn
"""

import argparse
import sys
import random
import turtle
import math
from Bio import SeqIO

__version__ = '2.0'

def get_argument():
    # Parsers
    parser = argparse.ArgumentParser(description = 'geneturtle', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_group_1 = parser.add_argument_group('Input and Output')
    parser_group_2 = parser.add_argument_group('Parameters')

    # Input and output
    parser_group_1.add_argument('-i', '--input', required = True, type = str, 
                                help = 'The file include the detail to draw the figure')
    parser_group_1.add_argument('-t', '--type', required = True, type = str, 
                                help = 'The file type of your input file, "gbk" and "list" are available')
    parser_group_1.add_argument('-o', '--output', required = False, type = str, default = 'gene_turtle',
                              help = 'Output file')
    
    # Parameters
    parser_group_2.add_argument('-c', '--color_mode', required = False, type = str, default = 'preset',
                                help = 'We have three mode for color assignment for every gene group, you can choose "random", "preset", or upload the color file set by yourself. [default: preset]')
    parser_group_2.add_argument('-f', '--force', action = 'store_true', help = 'In case the sequence is too long and you still want to run it.')
    parser_group_2.add_argument('-v', '--version', action = 'version', version = 'geneturtle v' + __version__, 
                        help = 'Show version number and exit')
    return parser

def feature_extraction(inputfile, seq_type):
    # translate the input file to a table contain 5 values, gene name, start, end, strand, and group
    features = []
    if seq_type == 'list':
        # read input list
        feature_file = open(inputfile, 'rt').readlines()
        # iterate detail of every single gene
        for gene in feature_file:
            gene = gene.strip().split(' ')
            # find the start, end, strand, and in case the group contain more than one word
            if type(gene[1]) == str:
                gene[1] = eval(gene[1])
                gene[2] = eval(gene[2])
            gene[4] = ' '.join(gene[4:])
            gene = gene[:5]
            features.append(gene)
    elif seq_type == 'gbk':
        # read the input genbank file
        for record in SeqIO.parse(inputfile, "gb"):
            # the gene information is within features object
            for gene in record.features:
                # only extract the information of every cds
                if gene.type == 'CDS':
                    # in case some gene do not have a name, use locus tag instead.
                    try:
                        name = gene.qualifiers['gene'][0]
                    except:
                        name = gene.qualifiers['locus_tag'][0]
                    # extract the start, end, strand and product
                    start = gene.location.start.numerator
                    end = gene.location.end.numerator
                    if gene.location.strand == 1:
                        strand = '+'
                    else:
                        strand = '-'
                    try:
                        product = gene.qualifiers['product'][0]
                    except:
                        product = 'hypothetical protein'
                    features.append([name, start, end, strand, product])
    # in case some inputfile not start by 1
    if features[0][1] != 1:
        normalize = features[0][1] - 1
        for feature in features:
            feature[1] -= normalize
            feature[2] -= normalize
    return features

def short_cut(features):
    size = features[-1][2] // 10000 + 1
    # how many time you want to make it short,this number could be adjust by yourself here
    simplify_time = size / 4
    for feature in features:
        feature[1] = feature[1] / simplify_time
        feature[2] = feature[2] / simplify_time
    return features, simplify_time
    

def cluster_the_annotation(features, seq_type, colormode):
    # clustal all the features by given colormode
    clusters_with_color = {}
    if colormode == 'preset':
        if seq_type == 'gbk':
            features = concise_gbk_feature(features)
        #the default color group is more soft
        colors = [[231, 219, 202], [184, 241, 204], [184, 241, 237], [241, 241, 184], [241, 204, 184], [217, 184, 241], [241, 184, 241], [255, 155, 106], [221, 255, 149], [184, 211, 143]]
        #provide another group of color, which is more bright
        #colors = [[102, 204, 204], [204, 255, 102], [255, 153, 204], [255, 0, 51], [255, 153, 102], [204, 51, 153], [255, 255, 102], [153, 204, 51], [102, 102, 153], [102, 102, 51]]
        #collect all groups name
        clusters = []
        for feature in features:
            if feature[4] not in clusters:
                clusters.append(feature[4])
            else:
                continue
        # only ten preset color, so this is a threshold
        if len(clusters) > 10:
            print('Too many functional annotation cluster, please consice it or choose other colormode')
            sys.exit(0)
        # random choose the color
        choosed_color = []
        for i in range(len(clusters)):
            while True:
                color = random.choice(colors)
                if color not in choosed_color:
                    choosed_color.append(color)
                    break
            # assign selected colors to every group
            clusters_with_color[clusters[i]] = color
    elif colormode == 'random':
        if seq_type == 'gbk':
            features = concise_gbk_feature(features)
        clusters = []
        for feature in features:
            if feature[4] not in clusters:
                clusters.append(feature[4])
            else:
                continue
        choosed_color = []
        for i in range(len(clusters)):
            while True:
                color = [random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)]
                if color not in choosed_color:
                    choosed_color.append(color)
                    break
            clusters_with_color[clusters[i]] = color
    else:
        # the colormode[1] should be the extra group and color sheet, for example: Transposase 0,255,255.
        color_file = open(colormode, 'rt').readlines()
        # assign the colors to every group
        for line in color_file:
            line = line.strip().split(' ')
            group_name, color = ' '.join(line[0:-1]), line[-1]
            clusters_with_color[group_name] = color.split(',')
    return clusters_with_color

def concise_gbk_feature(features):
    # this table need in time update
    group_preset = {
        'Transposase' : ['transposase'], 
        'Transcriptional regulator' : ['regulator'], 
        'Efflux transporter' : ['efflux transporter'], 
        'Hydrolase' : ['hydrolase'], 
        'Reductase' : ['reductase'], 
        'ATPase' : ['ATPase'], 
        'Recombinase' : ['recombinase'], 
        'Hypothetical protein' : ['hypothetical protein'], 
        'Mutase' : ['mutase'], 
        }
    # the amr pending is depend on the gene name, and others are product
    amr_preset = ['sul', 'tet', 'bla', 'cat']
    # pending the group
    for feature in features:
        # pending if the gene is a amr
        for amr in amr_preset:
            if amr in feature[0]:
                feature[4] = 'Antibiotic resistance'
        if feature[4] == 'Antibiotic resistance':
            continue
        else:
            # pending if a gene could be assigned to a group, else it will be noted as hp.
            positive_mark = 0
            for group_name, keywords in group_preset.items():
                for keyword in keywords:
                    if keyword in feature[4]:
                        feature[4] = group_name
                        positive_mark = 1
            if positive_mark == 0:
                feature[4] = 'Hypothetical protein'
    return features

def initiate(features, simplify_time):
    # count the total length
    total_length = features[-1][2]
    # set the window
    turtle.setup(1200, 800)
    # make it quick
    turtle.speed(0)
    # choose the color mode
    turtle.colormode(255)
    turtle.penup()
    # count the line number
    n_line = (total_length // 10000 + 1) / simplify_time
    # set the position of y axis
    y_axis = 0
    # we cannot draw a sequence longer than 40k,
    if n_line > 4:
        print('Your sequence is longer than 40, 000 bp, and there maybe some unpredictable errors, if you still want to do it, please add "-f" or "--force"')
        sys.exit(0)
    elif n_line == 4:
        y_axis = 300
    elif n_line == 3:
        y_axis = 250
    elif n_line == 2:
        y_axis = 200
    else:
        y_axis = 150
    # go to the start point
    turtle.goto(-520, y_axis)
    # set the pen size
    turtle.pensize(2)
    turtle.pendown()
    # we can give it a 20 start
    turtle.fd(20)
    return y_axis, n_line

def draw_line_in_edge(distance_to_move, length_for_one_line, y_axis, n_line):
    # This funcion is only aim to draw the line between the arrows
    # only draw the the line to 1000
    this_line = 1000 - length_for_one_line + distance_to_move
    # next line will draw the last part
    next_line = distance_to_move - this_line
    # draw the line and the 20 tail
    turtle.fd(this_line)
    turtle.fd(20)
    turtle.penup()
    # move the turtle to next line, pass this if only one line left
    if n_line >= 2:
        y_axis -= 150
        n_line -= 1
    turtle.goto(-520, y_axis)
    # draw the 20 head and rest part
    turtle.pendown()
    turtle.fd(20)
    turtle.fd(next_line)
    # the next line will return as the initial length of next line
    return y_axis, n_line, next_line

def draw_edge_fw(length, color, y_axis, n_line, length_for_one_line, gene_name):
    # only draw the arrow in this line to 1000
    this_line = 1000 - length_for_one_line + length
    next_line = length - this_line
    # draw the rectangle will no right barrier, fill in by selected color, the height of the rectangle is 40
    turtle.penup()
    turtle.fd(this_line)
    turtle.left(90)
    turtle.fd(20)
    turtle.left(90)
    turtle.pendown()
    turtle.fillcolor(eval(color[0]), eval(color[1]), eval(color[2]))
    turtle.begin_fill()
    turtle.fd(this_line)
    turtle.left(90)
    turtle.fd(40)
    turtle.left(90)
    turtle.fd(this_line)
    turtle.end_fill()
    turtle.left(90)
    turtle.penup()
    turtle.fd(20)
    turtle.right(90)
    # pending were is the proper place to write the gene name.
    if this_line >= next_line:
        # if the gene is too short, then the name will not be write.
        if not this_line < 35:
            add_cite(gene_name, this_line)
            turtle.penup()
    # go to next line
    if n_line >= 2:
        y_axis -= 150
        n_line -= 1
    # here is -500 because we don't need to draw the 20 head
    turtle.goto(-500, y_axis)
    # start to draw the rest of the arrow
    turtle.left(90)
    turtle.fd(20)
    turtle.right(90)
    turtle.pendown()
    turtle.begin_fill()
    turtle.fd(next_line - 35)
    turtle.left(90)
    turtle.fd(15)
    turtle.right(135)
    turtle.fd(35/math.cos(math.pi/4))
    turtle.right(90)
    turtle.fd(35/math.cos(math.pi/4))
    turtle.right(135)
    turtle.fd(15)
    turtle.left(90)
    turtle.fd(next_line - 35)
    turtle.right(90)
    turtle.end_fill()
    # move the turtle to the loc of turtle
    turtle.penup()
    turtle.fd(20)
    turtle.right(90)
    turtle.fd(next_line)
    turtle.pendown()
    # pending were is the proper place to write the gene name.
    if this_line < next_line:
        # if the gene is too short, then the name will not be write.
        if not next_line < 35:
            add_cite(gene_name, next_line)
    return y_axis, n_line
        
def draw_edge_rev(length, color, y_axis, n_line, length_for_one_line, gene_name):
    # only draw the arrow in this line to 1000
    this_line = 1000 - length_for_one_line + length
    next_line = length - this_line
    # you have to make the triangle complete, and 35 for the triangle, 2 for a tiny tail to make it more beautiful.
    if this_line <= 35:
        next_line -= 37 - this_line
        this_line = 37
    # move to the end and draw the arrow reverse, to make sure there are no right barrier.
    turtle.penup()
    turtle.fd(this_line)
    turtle.left(90)
    turtle.fd(20)
    turtle.left(90)
    turtle.pendown()
    turtle.fillcolor(eval(color[0]), eval(color[1]), eval(color[2]))
    turtle.begin_fill()
    turtle.fd(this_line - 35)
    turtle.right(90)
    turtle.fd(15)
    turtle.left(135)
    turtle.fd(35/math.cos(math.pi/4))
    turtle.left(90)
    turtle.fd(35/math.cos(math.pi/4))
    turtle.left(135)
    turtle.fd(15)
    turtle.right(90)
    turtle.fd(this_line - 35)
    turtle.left(90)
    turtle.end_fill()
    turtle.penup()
    turtle.fd(20)
    turtle.right(90)
    # pending were is the proper place to write the gene name.
    if this_line >= next_line:
        # if the gene is too short, then the name will not be write.
        if not this_line < 35:
            add_cite(gene_name, this_line)
            turtle.penup()
    # go to next line
    if n_line >= 2:
        y_axis -= 150
        n_line -= 1
    # here is -500 because we don't need to draw the 20 head
    turtle.goto(-500, y_axis)
    # start to draw the rest of the arrow
    turtle.left(90)
    turtle.fd(20)
    turtle.right(90)
    turtle.pendown()
    turtle.begin_fill()
    turtle.fd(next_line)
    turtle.right(90)
    turtle.fd(40)
    turtle.right(90)
    turtle.fd(next_line)
    turtle.right(90)
    turtle.end_fill()
    turtle.penup()
    turtle.fd(20)
    turtle.right(90)
    turtle.fd(next_line)
    turtle.pendown()
    # pending were is the proper place to write the gene name.
    if this_line < next_line:
        # if the gene is too short, then the name will not be write.
        if not next_line < 35:
            add_cite(gene_name, next_line)
    return y_axis, n_line

def draw_fw_arrow(length, color):
    turtle.fillcolor(eval(color[0]), eval(color[1]), eval(color[2]))
    turtle.begin_fill()
    turtle.left(90)
    turtle.fd(20)
    turtle.right(90)
    # if the length of the gene is less than twice length of the right triangle head, then half of the length will be the rectangle, and the rest half is a triangle.
    if length >= 70:
        turtle.fd(length - 35)
    else:
        turtle.fd(length/2)
    turtle.left(90)
    turtle.fd(15)
    if length >= 70:
        turtle.right(135)
        turtle.fd(35/math.cos(math.pi/4))
        turtle.right(90)
        turtle.fd(35/math.cos(math.pi/4))
        turtle.right(135)
        turtle.fd(15)
        turtle.left(90)
    else:
        # calculate the angle and the length of the triangle
        radian = math.atan(length/2/35)
        angle = radian / math.pi * 180
        turtle.right(180 - angle)
        turtle.fd(35/math.cos(radian))
        turtle.right(2 * angle)
        turtle.fd(35/math.cos(radian))
        turtle.right(180 - angle)
        turtle.fd(15)
        turtle.left(90)
    if length >= 70:
        turtle.fd(length - 35)
    else:
        turtle.fd(length/2)
    turtle.right(90)
    turtle.fd(20)
    turtle.right(90)
    turtle.end_fill()
    turtle.penup()
    turtle.fd(length)
    turtle.pendown()
    
def draw_rev_arrow(length, color):
    turtle.fillcolor(eval(color[0]), eval(color[1]), eval(color[2]))
    turtle.begin_fill()
    # if the length of the gene is less than twice length of the right triangle head, then half of the length will be the rectangle, and the rest half is a triangle.
    if length >= 70:
        turtle.left(45)
        turtle.fd(35/math.cos(math.pi/4))
        turtle.right(135)
    else:
        radian = math.atan(length/2/35)
        angle = radian / math.pi * 180
        turtle.left(90 - angle)
        turtle.fd(35/math.cos(radian))
        turtle.right(180 - angle)
    turtle.fd(15)
    turtle.left(90)
    if length >= 70:
        turtle.fd(length - 35)
    else:
        turtle.fd(length/2)
    turtle.right(90)
    turtle.fd(40)
    turtle.right(90)
    if length >= 70:
        turtle.fd(length - 35)
    else:
        turtle.fd(length/2)
    turtle.left(90)
    turtle.fd(15)
    if length >=70:
        turtle.right(135)
        turtle.fd(35/math.cos(math.pi/4))
        turtle.right(135)
    else:
        turtle.right(180 - angle)
        turtle.fd(35/math.cos(radian))
        turtle.right(90 + angle)
    turtle.end_fill()
    turtle.penup()
    turtle.fd(length)
    turtle.pendown()
    
def add_cite(gene_name, length):
    loc_to_write = length
    # trying to put the name in the center
    if length / len(gene_name) > 20:
        loc_to_write = length/4*3
    # write the gene name 30 inch under the arrow
    turtle.penup()
    turtle.bk(loc_to_write)
    turtle.left(90)
    turtle.bk(50)
    turtle.write(gene_name, font = ('Times new roman', 10, 'normal'))
    # back to the origin place
    turtle.fd(50)
    turtle.right(90)
    turtle.fd(loc_to_write)
    turtle.pendown()

def figure_legends(cluster_the_annotation, length_for_one_line, simplify_time):
    # move the turtle to the bottom of the paper, left side
    turtle.penup()
    turtle.right(90)
    turtle.fd(100)
    turtle.left(90)
    turtle.bk(length_for_one_line)
    turtle.pendown()
    # draw boxes for every group, 5 per line max.
    count = 0
    for cluster, color in cluster_the_annotation.items():
        count += 1
        # 5 annotation for each line
        if count == 6:
            turtle.penup()
            turtle.bk(650)
            turtle.right(90)
            turtle.fd(50)
            turtle.left(90)
            turtle.pendown()
        draw_a_box(cluster, color)
    # draw the scale
    turtle.penup()
    turtle.goto(450, -300)
    turtle.pendown()
    turtle.left(90)
    turtle.fd(5)
    turtle.bk(5)
    turtle.right(90)
    if 1 <= simplify_time <= 2:
        turtle.fd(100 / simplify_time)
        turtle.left(90)
        turtle.fd(5)
        turtle.bk(5)
        turtle.right(90)
        turtle.penup()
        turtle.bk(55 / simplify_time)
        turtle.right(90)
        turtle.fd(20)
        turtle.write('1 kb', font = ('Times new roman', 10, 'normal'))
    elif 2 < simplify_time <= 4:
        turtle.fd(200 / simplify_time)
        turtle.left(90)
        turtle.fd(5)
        turtle.bk(5)
        turtle.right(90)
        turtle.penup()
        turtle.bk(110 / simplify_time)
        turtle.right(90)
        turtle.fd(20)
        turtle.write('2 kb', font = ('Times new roman', 10, 'normal'))
    else:
        print('the length of the sequence is not suitable for geneturtle.')
        sys.exit(0)

def draw_a_box(cluster, color):
    turtle.fillcolor(eval(color[0]), eval(color[1]), eval(color[2]))
    turtle.begin_fill()
    turtle.fd(25)
    turtle.left(90)
    turtle.fd(25)
    turtle.left(90)
    turtle.fd(25)
    turtle.left(90)
    turtle.fd(25)
    turtle.end_fill()
    turtle.penup()
    turtle.left(90)
    turtle.fd(30)
    turtle.write(cluster, font = ('Times new roman', 12, 'normal'))
    turtle.fd(100)
    turtle.pendown()

def main():
    print('If you have any questions or suggestions for geneturtle, please contact Genglin Guo, e-mail: 2019207025@njau.edu.cn')
    # parse the input
    args = get_argument().parse_args()
    # translate the input file to a table contain 5 values, gene name, start, end, strand, and group
    features = feature_extraction(args.input, args.type)
    # make it short
    simplify_time = 1
    if args.force:
        features, simplify_time = short_cut(features)
    # clustal all the features by given colormode
    clusters_with_color = cluster_the_annotation(features, args.type, args.color_mode)
    # set the initiate
    y_axis, n_line = initiate(features, simplify_time)
    # count the length for one line, if it > 1000, change to next line.
    length_for_one_line = 0
    # trace the location of the turtle
    the_loc_of_turtle = 1
    # draw the arrow turtle by turtle
    for feature in features:
        if feature[1] > the_loc_of_turtle:
            # 10bp (multiply by simplify_time) for every inch
            distance_to_move = (feature[1] - the_loc_of_turtle) / 10
            length_for_one_line += distance_to_move
            # it the length for one line was over 1000 after draw the line
            if length_for_one_line >= 1000:
               y_axis, n_line, length_for_one_line = draw_line_in_edge(distance_to_move, length_for_one_line, y_axis, n_line)
            else:
                # if not over 1000, just draw the line
                turtle.fd(distance_to_move)
        # sometimes the genes will have some overlap
        elif feature[1] < the_loc_of_turtle:
            # let the turtle go back
            distance_to_move = the_loc_of_turtle - feature[1]
            turtle.penup()
            turtle.bk(distance_to_move)
            # make the turtle back to the origin status
            turtle.pendown()
            length_for_one_line -= distance_to_move
        # no matter what, we we have drawed one arrow, the turtle will located at the end of the arrow
        the_loc_of_turtle = feature[2]
        # find the correlated color
        color = clusters_with_color[str(feature[4])]
        # re-new the length for one line
        length = ((feature[2] - feature[1]) + 1) / 10
        length_for_one_line += length
        # 35 will draw a right triangle, left the right triangle to next line is not beautiful.
        if length_for_one_line > 1035:
            if feature[3] == '+':
                y_axis, n_line = draw_edge_fw(length, color, y_axis, n_line, length_for_one_line, feature[0])
            else:
                y_axis, n_line = draw_edge_rev(length, color, y_axis, n_line, length_for_one_line, feature[0])
            # reset the length for one line and minus the length we already draw
            length_for_one_line = length_for_one_line - 1000
        # if the length will not exceed to much, the arrow will still draw in this line, and move the turtle to next line.
        elif 1000 <= length_for_one_line <= 1035:
            if feature[3] == '+':
                draw_fw_arrow(length, color)
            else:
                draw_rev_arrow(length, color)
            # 20 inch tiny tail
            turtle.fd(20)
            # move to next line
            turtle.penup()
            if n_line >= 2:
                y_axis -= 150
                n_line -= 1
            turtle.goto(-520, y_axis)
            turtle.pendown()
            # 20 inch tiny head
            turtle.fd(20)
            # reset the length for one line
            length_for_one_line = 0
        else:
            # this is the normal situation, simply draw a arrow
            if feature[3] == '+':
                draw_fw_arrow(length, color)
            else:
                draw_rev_arrow(length, color)
            if not length < 35:
                add_cite(feature[0], length)
    turtle.fd(20)
    figure_legends(clusters_with_color, length_for_one_line, simplify_time)
    turtle.hideturtle()
    ts = turtle.getscreen()
    file_name = args.output + '.eps'
    ts.getcanvas().postscript(file = file_name)
    
main()
