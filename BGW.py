#!/usr/bin/env python3

import random
import math
import statistics
from scipy.special import binom
from tkinter import *
from itertools import accumulate

def poisson_p(rate ,k):
    return (rate ** k) * math.exp(-rate)/math.factorial(k)

def poisson_phi(rate, t):
    return math.exp(rate * (t-1))

def poisson_phi_inv(rate, x):
    if x == 0:
        return float("-inf")
    return (rate + math.log(x))/rate

def binomial_p(d, p, k):
    return binom(d, k) * (p ** k) * ((1-p) ** (d-k))

def binomial_phi(d, p, t):
    return ((1-p)+p*t) ** d

def binomial_phi_inv(d, p, x):
    return (math.pow(x, 1/d) + p-1) / p

def geometric_p(p, k):
    return ((1-p) ** (k-1)) * p

def geometric_phi(p, t):
    return (1-p) / (1-p*s)

def geometric_phi_inv(p, x):
    return (x+p-1)/(p*x)

def smallest_fixed_point(f, err):
    x = 0.
    while abs(f(x)-x) > err and x < 1:
        x += err
    return x

def smallest_fixed_point2(phi_inv, err):
    x = .99
    while abs(phi_inv(x)-x) > err:
        x = phi_inv(x)
    return x

def estimate_mean(phi, dt):
    return (phi(1)-phi(1-dt))/dt

def p(k):
    return binomial_p(3,.5, k)

def phi(t):
    return binomial_phi(3,.5,t)

def phi_inv(x):
    return binomial_phi_inv(3,.5,x)

class CDF:
    def __init__(self, p):
        self.p = p
        self.vals = list(accumulate((self.p(i) for i in range(20))))

    def make_more_vals(self):
        n = len(self.vals)
        inc = lambda x: x + self.vals[-1]
        new_vals = map(inc, accumulate((self.p(i) for i in range(n, 2*n))))
        self.vals.extend(new_vals)

    def inv(self, u):
        while u > self.vals[-1]:
            self.make_more_vals()
        imax = len(self.vals)-1
        imin = 0
        while imin < imax:
            i = int(imin + (imax-imin)/2)
            if u > self.vals[i]:
                imin = i+1
            else:
                imax = i
        return imin

    def inv2(self, u):
        for i in range(len(self.vals)):
            if u < self.vals[i]:
                return i

    def random(self):
        return self.inv(random.random())

class BGWTree:
    class Node:
        def __init__(self, parent):
            self.children = []
            self.parent = parent
    
        def make_children(self, cdf):
            n = cdf.random()
            self.children = [BGWTree.Node(self) for x in range(n)]

        def num_children(self):
            return len(self.children)

    def __init__(self, max_depth, p):
        self.root = BGWTree.Node(None)
        self.generations = [[self.root]]
        self.p = p
        self.cdf = CDF(p)
        self.max_depth = max_depth
        generation_number = 0
        while generation_number < self.max_depth:
            if not self.generations[generation_number]:
                # extinction
                break
            next_generation = []
            for node in self.generations[generation_number]:
                node.make_children(self.cdf)
                next_generation.extend(node.children)
            self.generations.append(next_generation)
            generation_number += 1

    def generation_size(self, n):
        if 0 <= n and n < self.num_generations():
            return len(self.generations[n])
        return 0

    def num_generations(self):
        return len(self.generations)

    def goes_extinct(self):
        return len(self.generations)-1 < self.max_depth

    def has_d_ary_subtree(self, d):
        return self.has_d_ary_subtree_helper(d, self.root, self.max_depth)

    def has_d_ary_subtree_helper(self, d, node, height):
        if height==0:
            return True
        if node.num_children() < d:
            return False
        n = 0
        for child in node.children:
            if self.has_d_ary_subtree_helper(d, child, height-1):
                n += 1
            if n == d:
                return True
        return False

class BGWDemo(Frame):
    def __init__(self, master=None, p=None, phi=None, phi_inv=None, max_depth=7):
        Frame.__init__(self, master)
        self.p = p
        self.max_depth = max_depth
        self.num_sampled = 0
        self.num_extinct = 0
        self.num_with_binary_subtree = 0
        self.generation_sizes = [[] for i in range(max_depth+1)]
        self.generation_sizes_conditioned_extinct = [[] for i in range(max_depth+1)]
        self.extinct_phi = None
        self.given_extinction_extinct_p = None
        self.mean = -1
        self.mean_given_extinction = -1
        self.extinct_p = None
        if phi:
            self.extinct_p = round(smallest_fixed_point(phi,.0001), 3)
            self.mean = estimate_mean(phi, .0001)
            self.theoretical_generation_sizes = [self.mean ** i for i in range(max_depth+1)]
        if phi_inv:
            self.extinct_p = round(smallest_fixed_point2(phi_inv,.0001), 3)
        if phi and self.extinct_p > 0:
            self.extinct_phi = lambda x: phi(self.extinct_p * x)/ self.extinct_p
            self.given_extinction_extinct_p = round(smallest_fixed_point(self.extinct_phi,.0001), 3)
            self.mean_given_extinction = estimate_mean(self.extinct_phi, .0001)
            self.theoretical_generation_sizes_conditioned_extinct = [self.mean_given_extinction ** i for i in range(max_depth+1)]

        self.canvas_frame = Frame(master)
        self.buttons_frame = Frame(master)
        self.canvas = Canvas(self.canvas_frame, highlightthickness=0)
        self.master.bind("<Configure>", self.on_resize)
        self.new_tree_button = Button(self.buttons_frame, text="Draw",
                command=self.new_tree_and_redraw) 
        self.quit_button = Button(self.buttons_frame, text="Quit", command=root.destroy)
        self.extinct_label = Label(self.buttons_frame, text="Goes extinct:")
        self.has_binary_subtree_label = Label(self.buttons_frame, text="Has binary subtree:")
        self.percent_with_binary_subtree_label = Label(self.buttons_frame,
                text="Percent with binary subtree:")
        self.extinct_percent_label = Label(self.buttons_frame, text="Percent extinct:")
        self.extinct_theoretical_label = Label(self.buttons_frame, text="Theretical prob of extinction:")
        self.num_sampled_label = Label(self.buttons_frame, text="Sampled:")
        self.draw_ten_button = Button(self.buttons_frame, text="Draw 10",
                command=self.sample(10))
        self.draw_hundred_button = Button(self.buttons_frame, text="Draw 100",
                command=self.sample(100))
        self.draw_thousand_button = Button(self.buttons_frame, text="Draw 1000",
                command=self.sample(1000))
        self.draw_ten_thousand_button = Button(self.buttons_frame, text="Draw 10000",
                command=self.sample(10000))
        self.generation_sizes_label = Label(self.buttons_frame, text="Avg generation sizes:") 
        self.generation_sizes_extinct_label = Label(self.buttons_frame, text="Extinct avg generation sizes:") 
        self.theoretical_generation_sizes_label = Label(self.buttons_frame, text="Theoretical generation size means:")
        self.mean_label = Label(self.buttons_frame, text="Theoretical mean:")
        self.mean_given_extinction_label = Label(self.buttons_frame, text="Theoretical mean given extinction:")
        self.theoretical_generation_sizes_given_extinction_label = Label(self.buttons_frame, text="Theoretical generation sizes given extinction:")
        
        self.num_sampled_label.pack(side=TOP)
        self.mean_label.pack(side=TOP)
        self.mean_given_extinction_label.pack(side=TOP)
        self.generation_sizes_label.pack(side=TOP)
        self.theoretical_generation_sizes_label.pack(side=TOP)
        self.generation_sizes_extinct_label.pack(side=TOP)
        self.theoretical_generation_sizes_given_extinction_label.pack(side=TOP)
        self.extinct_label.pack(side=TOP)
        self.extinct_percent_label.pack(side=TOP)
        self.extinct_theoretical_label.pack(side=TOP)
        self.has_binary_subtree_label.pack(side=TOP)
        self.percent_with_binary_subtree_label.pack(side=TOP)
        self.new_tree_button.pack(side=LEFT)
        self.draw_ten_button.pack(side=LEFT)
        self.draw_hundred_button.pack(side=LEFT)
        self.draw_thousand_button.pack(side=LEFT)
        self.draw_ten_thousand_button.pack(side=LEFT)
        self.quit_button.pack(side=LEFT)
        self.canvas.pack(expand=True, fill=BOTH)
        self.canvas_frame.pack(side=TOP, fill=BOTH, expand=True)
        self.buttons_frame.pack(side=BOTTOM)
        self.new_tree()

    def on_resize(self, event):
        self.draw_tree()

    def new_tree(self):
        self.num_sampled += 1
        self.tree = BGWTree(p=self.p, max_depth=self.max_depth)
        if self.tree.goes_extinct():
            self.num_extinct += 1
        if self.tree.has_d_ary_subtree(2):
            self.num_with_binary_subtree += 1
        for i in range(self.tree.max_depth+1):
            self.generation_sizes[i].append(self.tree.generation_size(i))
            if self.tree.goes_extinct():
                self.generation_sizes_conditioned_extinct[i].append(self.tree.generation_size(i))

    def new_tree_and_redraw(self):
        self.new_tree()
        self.draw_tree()

    def sample(self, n):
        def f():
            for i in range(n):
                self.new_tree()
            self.draw_tree()
        return f

    def draw_tree(self):
        if not self.num_sampled:
            return
        self.canvas.delete("all")
        tree = self.tree
        num_generations = tree.num_generations()
        canvas_width = int(self.canvas.winfo_width())
        canvas_height = int(self.canvas.winfo_height())
        tree.root.xmin = 0
        tree.root.xmax = canvas_width
        for i in range(0, num_generations):
            generation = tree.generations[i]
            for node in generation:
                num_children = node.num_children()
                node_width = node.xmax-node.xmin
                for j in range(num_children):
                    node.children[j].xmin = node.xmin + node_width/num_children*j 
                    node.children[j].xmax = node.xmin + node_width/num_children*(j+1)
                node.x = int(node.xmin + node_width/2)
                node.y = int((canvas_height/(num_generations+1))*(i+1))
                node.r = 3
                color = "black" if num_children > 0 or i == self.max_depth else "red"
                self.canvas.create_oval(node.x-node.r, node.y-node.r,
                        node.x+node.r, node.y+node.r, fill=color)
                if node.parent:
                    self.canvas.create_line(node.x, node.y, node.parent.x,
                            node.parent.y)
        self.update_labels()

    def update_labels(self):
        tree = self.tree
        self.num_sampled_label.configure(text="Sampled: {}".format(self.num_sampled))
        self.mean_label.configure(text="Theoretical mean: {}".format(round(self.mean,3)))
        self.extinct_label.configure(text="Goes extinct: YES" if
                tree.goes_extinct() else "Goes extinct: NO")
        self.extinct_percent_label.configure(text="Percent extinct: {}".format(
            round(self.num_extinct/self.num_sampled,3)))
        if self.extinct_p:
            self.extinct_theoretical_label.configure(text="Theoretical prob of extinction: {}".format(round(self.extinct_p,3)))
        self.percent_with_binary_subtree_label.configure(
            text="Percent with binary subtree: {}".format(
            round(self.num_with_binary_subtree/self.num_sampled,3)))
        self.has_binary_subtree_label.configure(text="Has binary subtree: {}".format(
            "YES" if self.tree.has_d_ary_subtree(2) else "NO"))
        self.generation_sizes_label.configure(text="Avg generation sizes: {}".format(
            str([round(statistics.mean(x),3) for x in self.generation_sizes])))
        self.theoretical_generation_sizes_label.configure(text="Theoretical avg generation sizes: {}".format(
            str([round(x,3) for x in self.theoretical_generation_sizes])))
        self.mean_given_extinction_label.configure(text="Theoretical mean given extinction: {}".format(round(self.mean_given_extinction,3)))
        self.generation_sizes_extinct_label.configure(text="Extinct avg generation sizes: {}".format(
            str([round(statistics.mean(x),3) for x in self.generation_sizes_conditioned_extinct if x])))
        self.theoretical_generation_sizes_given_extinction_label.configure(text="Theoretical exctinct avg generation sizes: {}".format(str([round(x,3) for x in self.theoretical_generation_sizes_conditioned_extinct])))



root = Tk()
demo = BGWDemo(master=root, p=p, phi=phi, max_depth=9)
demo.mainloop()
