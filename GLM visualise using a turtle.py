# Author Queubee.
# 7 december
# pasting all code from the jupyter notebook into python so we can use Turtle

from turtle import *
from random import choice
from itertools import cycle

# Globals
MAIN_DATA_FILE = "CvixLerC9.loc"
LEADING_LINES_TO_SKIP = 7
OUTPUT_FILENAME = "output RF matrix of " + MAIN_DATA_FILE + ".txt"

with open(MAIN_DATA_FILE, 'r') as f:
    for i in range(LEADING_LINES_TO_SKIP):
        # just read the skippable lines and dont do anything with them
        f.readline()
        
    identifier_dict = dict()
    
    for line in f:
        
        # at the end of the data, some more data is displayed, separated by a newline character
        # so we end it here
        if line == "\n":
            break
        
        # The marker identifier and number are contained in lines containing a semicolon
        if ";" in line:
            # isolate the identifier and the number
            l = line.split(' ')
            id_name = l[0]
            number = l[3].replace("\n", "")
            marker_id = f"{id_name}" # old : f"{id_name}_{number}"
            
            # to make sure a marker doesn't appear twice, throw an error if it does.
            if marker_id in identifier_dict.keys():
                raise ValueError
            
            # if no error occurs, create an empty list as value for this marker.
            else:
                identifier_dict[marker_id] = list()
        
        # if no semicolon is found, this line should consist of data for the current marker
        else:
            # first we have to parse the line to remove whitespaces
            # we'll check if the line complies to the rule "start with a space", at least.
            if line.startswith("  "):
                # now add all characters of this line to the list of the current marker
                identifier_dict[marker_id] += [char for char in line if char in {"a", "b", "-"}]    # code candy


# failed chi-square section
markers_to_reject = []

with open('chi_squares voor ' +MAIN_DATA_FILE+ ".txt", 'w') as out:

    out.write(f"marker\tchi-squared\tobserved_a\texpected_a\tobserved_b\texpected_b,total_genotypes\n")
    for k, values in identifier_dict.items():
        total_genotypes = (len(values) - values.count('-'))

        observed_a = values.count('a')
        observed_b = values.count('b')
        expected_a = total_genotypes / 2    # we expect a:b to be 1:3 where b is ab,bb,ba
        expected_b =  total_genotypes / 2
        chi_squared = ((((observed_a - expected_a)**2) / expected_a ) + \
                      (((observed_b - expected_b)**2) / expected_b ))
        print(k, observed_a, expected_a, observed_b,
              expected_b, total_genotypes, chi_squared)
        out.write(f"{k}\t{round(chi_squared,2)}\t{observed_a}\t{expected_a}\t{observed_b}\t{expected_b},{total_genotypes}\n")
        # note: total_genotypes gebaseerd op a's en b's, -'s niet meerekenen.

        if chi_squared > 3.841: # ignore this marker.
            markers_to_reject.append((k, chi_squared))
            # print(k,'will be ignored with a chi_squared of',chi_squared)
 

with open('chi_squares_verworpen.txt', 'w') as out:
    for marker, chi_squared in markers_to_reject:
        identifier_dict.pop(marker)
        out.write(f"{marker}\t{chi_squared}\n")
        print('removed', marker, 'from dataset with chi-square:',round(chi_squared,2),)


    

def calculateRF(ref_marker, compare_marker):
    
    # set a variable
    recombinant_count = 0
    unused_count = 0 # -'s
    # iterate over both lists simultaneously whilst comparing values
    # 
    for ref_genotype, comp_genotype in zip(ref_marker, compare_marker):
        
        ## UNCERTAIN IS IT Ref == COMP OR REF != COMP OR REF = B COMP = A??
        ## 
        if (ref_genotype == 'b' and comp_genotype == 'a') or (
            comp_genotype == 'b' and  ref_genotype == 'a'):
            recombinant_count += 1
        if ref_genotype == '-' or comp_genotype == '-':
            unused_count += 1
    
    return ((recombinant_count / (len(ref_marker) - unused_count)) * 100)


# iterate over each key, value pair 

recombinant_factor = dict()

recombinant_factors_try2 = dict()


# reference and compare_to markers.
for ref_marker, ref_genotype in identifier_dict.items():
    # and again, so we form a matrix 
    # NOTE; THIS MATRIX WONT HAVE THEIR AXES NICELY SORTED!!!!!!!!!!!
    recombinant_factors_try2[f"{ref_marker}"] = dict()
    for comp_marker, comp_genotype in identifier_dict.items():
        
        # skip same keys and skip previously already matched keys
        if (not ref_marker == comp_marker):
            recombinant_factor[f"{ref_marker}_VS_{comp_marker}"]  = calculateRF(ref_genotype, comp_genotype)

            recombinant_factors_try2[f"{ref_marker}"][f"{comp_marker}"] = calculateRF(ref_genotype, comp_genotype)
                                                  
for i, (k, v) in enumerate(recombinant_factor.items()):
    if i > 10:
        break
    
    print(k,v)
print(len(recombinant_factor))

for i, (k, v) in enumerate(recombinant_factors_try2.items()):
    if i > 10:
        break
    
    print(k,list(v.items())[:2])
print(len(recombinant_factors_try2))


# store to tsv

sorted_keys = list(recombinant_factors_try2.keys())
sorted_keys.sort()

factors = recombinant_factors_try2

with open (OUTPUT_FILENAME, 'w') as out:
    # set headers
    out.write("X\t") # topleft cell
    for column_header in sorted_keys:
        out.write(f"{column_header}\t")
    out.write("\n")
    
    # set leftmost cell key/marker-name
    for row_header in sorted_keys:
        out.write(f"{row_header}\t")
        for key in sorted_keys:
            if (not key == row_header):
                
                recombinant_distance = factors[row_header][key]
                out.write(f"{recombinant_distance}\t")
            else:
                # there's no distance between marker x with itself
                out.write("X\t")
        out.write("\n")

            


class Marker():
    def __init__(self, parent, dist_to_parent, marker_name):
        self.parent = parent
        self.distance = dist_to_parent
        self.name = marker_name
        
    def __lt__(self, other):
        return self.distance < other.distance

    def __str__(self):
        return f"{self.name}:\t\t{self.distance}"
    
class GeneLinkageMap():
    """ h"""
    def __init__(self, factor_dict):
        """ H"""
        self.markers = []
        self.factors = factor_dict
        self.madeMarkers = False
        self.sortedMarkers = False
        self.last = None
        self.now = None

    def setup(self):
        self.find_furthest_apart()
        self.create_markers()
    
        
    def create_markers(self):
        if self.madeMarkers:
            print('creating marker objects with ',self.primary_marker,
                  'as primary marker')
            self.markers.append(self.primary_marker)
        
            for marker_name, distance_to_parent in self.factors[self.primary_marker.name].items():
                marker = Marker(self.primary_marker.name, distance_to_parent, marker_name)
                self.markers.append(marker)
            
            print('sorting markers by distance to primary marker')
            self.markers.sort()
            self.sortedMarkers = True
        else:
            print("theres no markers to create")

        
        
    def display_map(self):
        for m in self.markers:
            print(m)

    def test_for_duplicates(self):
        print(len(self.markers), len(set(self.markers)))
        assert (len(self.markers) == len(set(self.markers)))
        print('no duplicates found')  
        
        
    def find_furthest_apart(self):
        print('finding 2 markers that are furthest apart')
        largest_distance = 0
        largest_1 = ''
        largest_2 = ''
        # the dictionary will look like this
        # { marker1 : {marker2 : RFcount, marker3: RFcount },
        #   marker2 : { marker1: RFcount, marker3: RFcount }}
        for marker in self.factors:
            for another_marker in self.factors:
                if (not another_marker == marker):
                    if self.factors[marker][another_marker] > largest_distance:
                       largest_distance = self.factors[marker][another_marker]
                       print(largest_distance)
                       largest_1 = marker
                       largest_2 = another_marker

        print('primary:',largest_1, '\nfurthest:',largest_2, '\ndistance:',
              largest_distance)
        if (not largest_1 == '' and not largest_2 == ''):
            self.primary_marker = Marker(largest_1, 0, largest_1)
            self.furthest_marker = Marker(largest_1, largest_distance, largest_2)
            
            self.madeMarkers = True
        else:
            print('no markers were found')
        
    def findNext(self):
        largest_distance = 0
        for marker in self.markers:
            pass

    def generate_mapchart(self):

        if not self.sortedMarkers:
            raise Exception( " Markers aren't sorted (yet) " )
        
        with open("MapchartValues - " + MAIN_DATA_FILE +" - .txt", 'w') as out:
            out.write("; a comment \n\n")
            out.write("some_group\n")

            for marker in self.markers:
                out.write(f"{marker.name}\t{marker.distance}\n")

            
    

    def display_turtlewise(self, coloring="RAINBOW", shapeshifter=True,
                           dumbsized=True):

        # todo alternate up down
        # todo tilt graph 90deg

        threshold = 1
        current = -100
        # # # # # # # # # # easter egg stuff # # # # # # # # # # # # # # # # #
        colors_of_the_rainbow = ['violet','red','orange','yellow','green','blue',
                                 'indigo']
        rainbow_colors = cycle(colors_of_the_rainbow)
        
        shapes = ["arrow", "turtle", "circle", "square", "triangle"]
        
        sizes = [size for size in range(3, 25)]

        c = next(rainbow_colors)
        for i, m in enumerate(self.markers):
            
            t = Turtle()

            # easter egg that colors the markers by the colors
            # of the rainbow
            if coloring == "RAINBOW":
                if round(i/len(self.markers)*100, 0) % len(
                               colors_of_the_rainbow) == 0:
                    c = next(rainbow_colors)
                # print(colors, c)
                t.color(c)
                
            # easter egg that colors markers randomly
            elif coloring == "RANDOM":
                t.color(choice(colors_of_the_rainbow))
                
            # easter egg that colors markers based on distance to parent
            elif coloring == "DISTANCE":
                fury_level =  (m.distance / 10)
                if fury_level > 0:
                    t.color('green')
                if fury_level > 2:
                    t.color('yellow')
                if fury_level > 3:
                    t.color('orange')
                if fury_level > 4:
                    t.color('red')
                if fury_level > 5:
                    t.color('brown')
                    
            # easter egg that changes shapes to flicker
            # randomly when a turtle is 'spawned'
            if shapeshifter:
                t.shape(choice(shapes))
                
            # easter egg that changes size to flicker
            # randomly when a turtle is spawned
            if dumbsized:
                t.turtlesize(choice(sizes))
                
            t.speed(0)
            t.hideturtle()
            t.left(90)
            t.forward(250)
            t.left(180)
            
            t.forward(9*m.distance)

            t.right(90)
            t.forward(80)

            t.goto(-150, (i*-12) + 250)
            t.up()
            t.goto(-290, (i*-12) + 245)
            t.down()
            t.write(m.name + " | " +f"{round(m.distance,2)}")

            current = m.distance
        ts = t.getscreen()
        ts.getcanvas().postscript(file="Gene Linkage Map with Python Turtle try 1.eps")
        
    

g = GeneLinkageMap(recombinant_factors_try2)
g.setup()
#g.display_map()
g.test_for_duplicates()
g.generate_mapchart()
g.display_turtlewise("RAINBOW", True, True)

