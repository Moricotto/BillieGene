from enum import Enum
from random import randint
from numpy.random import choice


from graphviz import Digraph
from PIL import Image

class Genes:
    R = 0
    r = 1

class Genotype(Enum):
    RR = 0
    Rr = 1
    rr = 2

class Phenotype(Enum):
    DOM = 0
    REC = 1

Genotype_to_Phenotype = {
    Genotype.RR: Phenotype.DOM,
    Genotype.Rr: Phenotype.DOM,
    Genotype.rr: Phenotype.REC
}

Genotype_to_Genes = {
    Genotype.RR: (Genes.R, Genes.R),
    Genotype.Rr: (Genes.R, Genes.r),
    Genotype.rr: (Genes.r, Genes.r)
}

Genes_to_Genotype = {
    (Genes.R, Genes.R): Genotype.RR,
    (Genes.R, Genes.r): Genotype.Rr,
    (Genes.r, Genes.R): Genotype.Rr,
    (Genes.r, Genes.r): Genotype.rr
}

class Node:
    def __init__(self, mate=None, children=[], required_phenotype=None, name=""):
        self.genotype = None
        self.required_phenotype = required_phenotype
        self.children = children
        self.mate = mate
        self.genotype_counts = [0, 0, 0]
        self.name = name


def mate(parent1, parent2):
    return Genes_to_Genotype[(Genotype_to_Genes[parent1.genotype][randint(0, 1)], Genotype_to_Genes[parent2.genotype][randint(0, 1)])]

PROBS = [0.25, 0.5, 0.25]  #RR, Rr, rr

class Tree:
    
    def __init__(self):
        self.nodes = []
        self.first_gen = []
        self.name = "A"
    
    def add_node(self, node, is_first_gen=False):
        self.nodes.append(node)
        #node.name = self.name
        #self.name = chr(ord(self.name) + 1)
        if is_first_gen:
            self.first_gen.append(node)

    def assign_genotypes(self):
        gen = self.first_gen
        next_gen = []
        while gen != []:
            for node in gen:
                if node.genotype is None: #no parents so not assigned
                    node.genotype = choice(list(Genotype), p=PROBS)
                    #print(f"Assigned {node.name} genotype {node.genotype}")
                if node.mate is not None and node.mate.genotype is None:
                    node.mate.genotype = choice(list(Genotype), p=PROBS)
                    #print(f"Assigned {node.mate.name} genotype {node.mate.genotype}")
                for child in node.children:
                    child.genotype = mate(node, node.mate)
                    child.phenotype = Genotype_to_Phenotype[child.genotype]
                    #print(f"Assigned {child.name} genotype {child.genotype} from parents {node.name} and {node.mate.name}")
                    next_gen.append(child)
            gen = next_gen
            next_gen = []
        #now we check for required phenotypes
        valid = True
        for node in self.nodes:
            if node.required_phenotype is not None:
                if Genotype_to_Phenotype[node.genotype] != node.required_phenotype:
                    #print(f"Node {node.name} has genotype {node.genotype} which does not match required phenotype {node.required_phenotype}")
                    valid = False
                    break
        if valid:
            for node in self.nodes:
                #print(f"Node {node.name} final genotype {node.genotype}")
                node.genotype_counts[node.genotype.value] += 1
                node.genotype = None  #reset for next attempt
        else:
            for node in self.nodes:
                node.genotype = None  #reset for next attempt

    def simulate(self, n):
        for _ in range(n):
            self.assign_genotypes()
    
    def convert_to_family_tree(self):
        family_tree = FamilyTree()
        for node in self.nodes:
            if node.mate is not None and family_tree.children.get(node.name) is None:
                family_tree.add_parents([child.name for child in node.children], node.name, node.mate.name)
            node.prob = (node.genotype_counts[Genotype.RR.value] + node.genotype_counts[Genotype.Rr.value]) / sum(node.genotype_counts)
            family_tree.set_trait_probability(node.name, node.prob)
        return family_tree


class FamilyTree:
    def __init__(self):
        self.parents = {}          # child -> set(parents)
        self.children = {}         # parent -> set(children)
        self.sibling_groups = []   # list of sibling lists
        self.trait_prob = {}       # person -> probability
        self.marriages = {}        # (parent1,parent2) -> marriage node

    # -----------------------------
    # Add parents for one or multiple children
    # -----------------------------
    def add_parents(self, children, parent1, parent2):
        if parent1 == parent2:
            raise ValueError("Child must have two different parents.")

        # Convert single child to list
        if isinstance(children, str):
            children = [children]

        for child in children:
            if child in self.parents and len(self.parents[child]) > 0:
                raise ValueError(f"Child '{child}' already has parents: {self.parents[child]}")
            self.parents[child] = {parent1, parent2}
            self.children.setdefault(parent1, set()).add(child)
            self.children.setdefault(parent2, set()).add(child)

        # Automatically add them as a sibling group
        if len(children) > 1:
            self.add_siblings(children)

    # -----------------------------
    # Add siblings manually (optional)
    # -----------------------------
    def add_siblings(self, siblings):
        if len(siblings) > 1:
            self.sibling_groups.append(siblings)

    # -----------------------------
    # Trait probabilities
    # -----------------------------
    def set_trait_probability(self, person, probability):
        if not (0 <= probability <= 1):
            raise ValueError("Probability must be between 0 and 1.")
        self.trait_prob[person] = probability

    def probability_to_color(self, p):
        if p is None:
            return "#ffffff"
        red = int(255 * p)
        green = 255 - int(80 * p)
        blue = 255 - int(80 * p)
        return f"#{red:02x}{green:02x}{blue:02x}"

    # -----------------------------
    # Generate Graph
    # -----------------------------
    def generate_graph(self, filename="family_tree"):
        g = Digraph(filename, format="png")
        g.attr(rankdir="TB", nodesep="0.5", ranksep="0.6")

        # Collect all people
        people = set(self.trait_prob.keys()) | set(self.parents.keys()) | set(self.children.keys())
        for kids in self.children.values():
            people.update(kids)
        for sibs in self.sibling_groups:
            people.update(sibs)

        # Add all person nodes
        for person in people:
            p = self.trait_prob.get(person)
            label = person if p is None else f"{person} ({p:.2f})"
            color = self.probability_to_color(p)
            g.node(person, label=label, shape="box", style="filled", fillcolor=color)

        # -----------------------------
        # Add marriage nodes and edges
        # -----------------------------
        # Group children by parents for single arrow
        parent_pairs = {}
        for child, parents in self.parents.items():
            pair = tuple(sorted(parents))
            parent_pairs.setdefault(pair, []).append(child)

        for (p1, p2), children in parent_pairs.items():
            marriage_node = f"marriage_{p1}_{p2}"
            self.marriages[(p1, p2)] = marriage_node
            g.node(marriage_node, shape="point", width="0.01", label="", style="invisible")
            g.edge(p1, marriage_node, dir="none")
            g.edge(p2, marriage_node, dir="none")

            # Connect first child from marriage node
            first_child = children[0]
            g.edge(marriage_node, first_child)

            # Connect other children via sibling lines (handled in sibling group below)

        # -----------------------------
        # Sibling lines (right->left)
        # -----------------------------
        for group in self.sibling_groups:
            with g.subgraph() as sg:
                sg.attr(rank="same")
                for i in range(len(group)-1):
                    left = group[i]
                    right = group[i+1]
                    sg.edge(f"{left}:e", f"{right}:w", dir="none")

        # -----------------------------
        # Render and display
        # -----------------------------
        outpath = g.render(filename, cleanup=True)
        print(f"[✔] Family tree saved as {outpath}")
        try:
            Image.open(outpath).show()
        except:
            pass


# -----------------------------
# Example usage
# -----------------------------
if __name__ == "__main__":
    tree = Tree()

    me = Node(required_phenotype=Phenotype.DOM, name="Me")
    mom = Node(name="Mom", children=[me])
    dad = Node(name="Dad", children=[me])
    mom.mate = dad
    dad.mate = mom
    cousin = Node(name="Cousin", required_phenotype=Phenotype.REC)
    aunt = Node(name="Aunt", children=[cousin])
    uncle = Node(name="Uncle", children=[cousin])
    aunt.mate = uncle
    uncle.mate = aunt
    grandma = Node(name="Grandma", children=[mom, aunt])
    grandpa = Node(name="Grandpa", children=[mom, aunt])
    grandma.mate = grandpa
    grandpa.mate = grandma
    tree.add_node(me)
    tree.add_node(mom)
    tree.add_node(dad)
    tree.add_node(aunt)
    tree.add_node(uncle)
    tree.add_node(grandma, is_first_gen=True)
    tree.add_node(grandpa)
    tree.simulate(100000)
    family_tree = tree.convert_to_family_tree()
    family_tree.generate_graph("family_tree_example")

