# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 15:54:04 2023

@author: Genglin Guo
"""

import sys
import turtle
import math
import random

def feature_extraction(inputfile):
    features = []
    feature_file = open(inputfile, 'rt')
    for line in feature_file:
        feature = []
        for i in line.split(' '):
            feature.append(i.strip())
        features.append(feature)
    if type(features[-1][2]) == str:
       for feature in features:
           feature[1] = eval(feature[1])
           feature[2] = eval(feature[2])
    if features[0][1] != 1:
        normalize = features[0][1] - 1
        for feature in features:
            feature[1] -= normalize
            feature[2] -= normalize
    return features

def initiate(features):
    if  type(features[-1][2]) == str:
        features[-1][2] = eval(features[-1][2])
        features[0][1] = eval(features[0][1])
    total_length = features[-1][2] - features[0][1]
    turtle.setup(1200, 800)
    turtle.colormode(255)
    turtle.penup()
    n_line = total_length // 10000 + 1
    y_axis = 0
    if n_line > 4:
        print('the input is too long to draw, gene_turtle is too weak to do this job, please try other tools')
        sys.exit(1)
    else:
        if n_line == 4:
            y_axis = 300
        elif n_line == 3:
            y_axis = 250
        elif n_line == 2:
            y_axis = 200
        else:
            y_axis = 150
    turtle.goto(-520, y_axis)
    turtle.pensize(2)
    turtle.pendown()
    turtle.fd(20)
    return y_axis, n_line

def cluster_the_annotation(features):
    clusters = []
    colors = [[231, 219, 202], [184, 241, 204], [184, 241, 237], [241, 241, 184], [241, 204, 184], [217, 184, 241], [241, 184, 241], [255, 155, 106], [221, 255, 149], [184, 211, 143]]
    for gene in features:
        if gene[4] not in clusters:
            clusters.append(gene[4])
        else:
            continue
    if len(clusters) > 10:
        print('Too many functional annotation cluster, please consice it')
        sys.exit(1)
    clusters_with_color = {}
    choosed_color = []
    for i in range(len(clusters)):
        while True:
            color = random.choice(colors)
            if color not in choosed_color:
                choosed_color.append(color)
                break
        clusters_with_color[clusters[i]] = color
    return clusters_with_color

def draw_line_in_edge(distance_to_move, length_for_one_line, y_axis, n_line):
    this_line = 1000 - length_for_one_line + distance_to_move
    next_line = distance_to_move - this_line
    turtle.fd(this_line)
    turtle.fd(20)
    turtle.penup()
    if n_line >= 2:
        y_axis -= 150
        n_line -= 1
    turtle.goto(-520, y_axis)
    turtle.pendown()
    turtle.fd(20)
    turtle.fd(next_line)
    return y_axis, n_line, next_line

def draw_fw_arrow(length, color):
    turtle.fillcolor(color[0], color[1], color[2])
    turtle.begin_fill()
    turtle.left(90)
    turtle.fd(20)
    turtle.right(90)
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
    turtle.fillcolor(color[0], color[1], color[2])
    turtle.begin_fill()
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
    if length / len(gene_name) > 20:
        loc_to_write = length/4*3
    turtle.penup()
    turtle.bk(loc_to_write)
    turtle.left(90)
    turtle.bk(50)
    turtle.write(gene_name, font = ('Times new roman', 10, 'normal'))
    turtle.fd(50)
    turtle.right(90)
    turtle.fd(loc_to_write)
    turtle.pendown()
        
def draw_edge_fw(length, color, y_axis, n_line, length_for_one_line, gene_name):
    this_line = 1000 - length_for_one_line + length
    next_line = length - this_line
    turtle.penup()
    turtle.fd(this_line)
    turtle.left(90)
    turtle.fd(20)
    turtle.left(90)
    turtle.pendown()
    turtle.fillcolor(color[0], color[1], color[2])
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
    if this_line >= next_line:
        if not this_line < 35:
            add_cite(gene_name, this_line)
    if n_line >= 2:
        y_axis -= 150
        n_line -= 1
    turtle.goto(-500, y_axis)
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
    turtle.penup()
    turtle.fd(20)
    turtle.right(90)
    turtle.fd(next_line)
    turtle.pendown()
    if this_line < next_line:
        if not next_line < 35:
            add_cite(gene_name, next_line)
    return y_axis, n_line
        
def draw_edge_rev(length, color, y_axis, n_line, length_for_one_line, gene_name):
    this_line = 1000 - length_for_one_line + length
    next_line = length - this_line
    if this_line <= 35:
        next_line -= 37 - this_line
        this_line = 37
    turtle.penup()
    turtle.fd(this_line)
    turtle.left(90)
    turtle.fd(20)
    turtle.left(90)
    turtle.pendown()
    turtle.fillcolor(color[0], color[1], color[2])
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
    if this_line >= next_line:
        if not this_line < 35:
            add_cite(gene_name, this_line)
    if n_line >= 2:
        y_axis -= 150
        n_line -= 1
    turtle.goto(-500, y_axis)
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
    if this_line < next_line:
        if not next_line < 35:
            add_cite(gene_name, next_line)
    return y_axis, n_line

def figure_legends(cluster_the_annotation, length_for_one_line):
    turtle.penup()
    turtle.right(90)
    turtle.fd(100)
    turtle.left(90)
    turtle.bk(length_for_one_line)
    turtle.pendown()
    count = 0
    for cluster, color in cluster_the_annotation.items():
        count += 1
        if count == 6:
            turtle.penup()
            turtle.bk(400)
            turtle.right(90)
            turtle.fd(50)
            turtle.left(90)
            turtle.pendown()
        draw_a_box(cluster, color)
    turtle.penup()
    turtle.goto(450, -300)
    turtle.pendown()
    turtle.left(90)
    turtle.fd(5)
    turtle.bk(5)
    turtle.right(90)
    turtle.fd(100)
    turtle.left(90)
    turtle.fd(5)
    turtle.bk(5)
    turtle.right(90)
    turtle.penup()
    turtle.bk(55)
    turtle.right(90)
    turtle.fd(20)
    turtle.write('1 kb', font = ('Times new roman', 10, 'normal'))
        
def draw_a_box(cluster, color):
    turtle.fillcolor(color[0], color[1], color[2])
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
    turtle.write(cluster, font = ('Times new roman', 10, 'normal'))
    turtle.fd(70)
    turtle.pendown()

def main():
    inputfile = sys.argv[-1]
    features = feature_extraction(inputfile)
    y_axis, n_line = initiate(features)
    clusters_with_color = cluster_the_annotation(features)
    length_for_one_line = 0
    the_loc_of_turtle = 1
    for gene in features:
        print(gene[0])
        if gene[1] > the_loc_of_turtle:
            distance_to_move = (gene[1] - the_loc_of_turtle) / 10
            length_for_one_line += distance_to_move
            if length_for_one_line >= 1000:
               y_axis, n_line, length_for_one_line = draw_line_in_edge(distance_to_move, length_for_one_line, y_axis, n_line)
            else:
                turtle.fd(distance_to_move)
        elif gene[1] < the_loc_of_turtle:
            distance_to_move = the_loc_of_turtle - gene[1]
            turtle.penup()
            turtle.bk(distance_to_move)
            turtle.pendown()
            length_for_one_line -= distance_to_move
        the_loc_of_turtle = gene[2]
        color = clusters_with_color[gene[4]]
        length = ((gene[2] - gene[1]) + 1) / 10
        length_for_one_line += length
        if length_for_one_line > 1035:
            if gene[3] == '+':
                y_axis, n_line = draw_edge_fw(length, color, y_axis, n_line, length_for_one_line, gene[0])
            else:
                y_axis, n_line = draw_edge_rev(length, color, y_axis, n_line, length_for_one_line, gene[0])
            length_for_one_line = length_for_one_line - 1000
        elif 1000 <= length_for_one_line <= 1035:
            if gene[3] == '+':
                draw_fw_arrow(length, color)
            else:
                draw_rev_arrow(length, color)
            turtle.fd(20)
            turtle.penup()
            if n_line >= 2:
                y_axis -= 150
                n_line -= 1
            turtle.goto(-520, y_axis)
            turtle.pendown()
            turtle.fd(20)
            length_for_one_line = 0
        else:
            if gene[3] == '+':
                draw_fw_arrow(length, color)
            else:
                draw_rev_arrow(length, color)
            if not length < 35:
                add_cite(gene[0], length)
    turtle.fd(20)
    figure_legends(clusters_with_color, length_for_one_line)
    turtle.hideturtle()
    ts = turtle.getscreen()
    ts.getcanvas().postscript(file="gene_map.eps")
    
main()

'''
基因box = [基因名，起始，终止，方向，功能注释]
例：[sul, 1, 1308, +, AMRG]

箭头以700bp为限
'''