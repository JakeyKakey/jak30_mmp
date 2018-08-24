from ete3 import Tree, TreeStyle, TextFace, NodeStyle
import pickle

ranked_lineage_filepath = '../bin/database/taxonomy/rankedlineage.dmp'

class TaxNode():

    def __init__(self, node_name, rank = None):
        self.node_name = node_name
        if rank == None:
            self.rank = ''
        else:
            self.rank = rank

        self.child_nodes = []
        self.parent_node = "ROOT"  # Should be overwritten in all instances but the root one.
        self.total_children = 0

    """"Adds a child node to the current node. """

    def addNode(self, node):
        self.child_nodes.append(node)

    """Check if the child node already exists """

    def childExists(self, node_name):
        for node in self.child_nodes:
            if node.node_name == node_name:
                return True
            else:
                return False

    """Gets specified child node based on its name """

    def getNode(self, node_name):
        for node in self.child_nodes:
            if node.node_name == node_name:
                return node

    """Sets the parent node of a child."""

    def setParent(self, parent_node):
        self.parent_node = parent_node

    """"Recursively counts the total amount of children nodes holds. Run once, after the tree is fully built. """

    def countChildren(self):
        for child in self.child_nodes:
            if(child.rank == 'Name'):
                self.total_children += 1

            self.total_children += child.countChildren()

        return self.total_children

    """Counts the most/least common match occurrences in a specific rank. Has two sub-methods,
        one for iterating over the tree and the other for printing out the collected results.
    
        Arguments:
        rank -- accepts 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Superkingdom'
        no_of_results -- number of results desired
        rootNode -- the absolute top node of the tree (shoddy coding I don't have time to fix)
        descending -- sort by most or least common. By default set to most common.
    """

    def calculateCommonnessByRank(self, rank, no_of_results, rootNode, descending = True):

        def iterate(node):
            for child in node.child_nodes:
                if(child.rank == rank):
                    percentage = (child.total_children / rootNode.total_children) * 100
                    list = []
                    list.append(percentage)
                    list.append(child.node_name)
                    results.append(list)

                iterate(child)

        def printResults():

            if descending:
                print("The most commonly occurring values for", rank, "are: ")
            else:
                print("The least commonly occurring values for", rank, "are: ")

            if len(results) > no_of_results:
                for i in range(no_of_results): # If this was any more complex than this, I would have refactored
                    result = results[i]        # to use a proper iterator.
                    print(str(result[0]),"% -",result[1])
            else:
                for result in results:
                    print(str(result[0]),"% -",result[1])

        no_of_results = int(no_of_results)

        results = []
        iterate(self)

        if descending:
            results = sorted(results, key=lambda x: x[0], reverse=True)
        else:
            results = sorted(results, key=lambda x: x[0])

        printResults()

"""Breaks up a rankedlineage.dmp line into a list of individual strings."""

def extract_data(lineage):
    lineage =' '.join(lineage.split())
    lineage = lineage.split("|")
    lineage.pop() # "" entry at the end
    return lineage

"""Assigns the ranks to the respective names in the lineage and fixes gaps, then adds to the tree."""

def prepare_for_ranked_tree(ranked_lineage, rootNode):
    ranks = ['Name', 'Species', 'Genus', 'Family', 'Order',
             'Class', 'Phylum', 'Kingdom', 'Superkingdom']

    combined = []

    count = 0
    for item in ranked_lineage: # A better approach would have been to use the double-iterate zip function
        list = []               # but I think this is readable enough.
        list.append(ranks[count])
        list.append(item)
        combined.append(list)
        count += 1

    combined = fixLineage(combined) # fixes the gaps in the lineage.
    build_ranked_tree(combined, rootNode)

"""Actual recursive generation of the tree."""

def build_ranked_tree(path, parent_node):
    while path:
        current_node = path.pop()
        current_node_name = current_node[1]
        current_node_rank = current_node[0]

        if (parent_node.getNode(current_node_name) and path):
            current_node = parent_node.getNode(current_node_name)
            build_ranked_tree(path, current_node)
        else:
            current_node = TaxNode(current_node_name,current_node_rank)
            parent_node.addNode(current_node)
            current_node.setParent(parent_node)
            build_ranked_tree(path, current_node)

"""Fixes the gaps and correctly matches the name to its respective rank """

def fixLineage(taxonomy_list):

    # So, this is doing a lot of obtuse looking direct array operations, but there's a reason for that.
    #
    # The list holds lists of two so taxonomy_list[0] will contain another list with
    # ["Candidatus Iainarchaeum andersonii","Name"], which should make some of what is happening more obvious.
    #
    # A ranked lineage may look like this:
    # Candidatus Iainarchaeum andersonii|	|Candidatus Iainarchaeum|	|	|	|Candidatus Diapherotrites|	|Archaea|
    #
    # This corresponds to Name|Species|Genus|Family|Order|Class|Phylum|Kingdom|Superkingdom
    # Species is empty which means Candidatus needs to be moved up there as it *is* the species.
    # This is so you can still have subspecies like Homo Neanderthalis > Homo Sapien, without a specific sub-species rank.
    #
    # There's data missing in Family, Order, Class, Kingdom which means this lineage just doesn't belong to a specific
    # subset and these need to be removed from the chain.

    name = taxonomy_list[0] # take out the Name
    remaining_lineage = taxonomy_list[1:9] # keep Species to Superkingdom here.

    last_empty = []
    to_remove = []

    # Moving up name.
    for item in remaining_lineage:
        if item[1] == ' ':
            last_empty = item
            to_remove.append(item)
        else:
            break

    # Removing all the gaps
    for item in remaining_lineage:
        if item[1] == ' ':
            to_remove.append(item)
    remaining_lineage = [x for x in remaining_lineage if x not in to_remove]

    #Associate the name with the last empty rank, and put that rank and name back onto the lineage.

    if len(last_empty) == 0:                # Strangely enough if not last_empty doesn't work, but this just puts
        remaining_lineage.insert(0, name)   # the name back at the beginning if there was nothing to move up
    else:
        last_empty[1] = name[1]
        remaining_lineage.insert(0, last_empty)
        remaining_lineage.insert(0, name)

    return remaining_lineage

"""Loads the rankedlineage.dmp taxonomy file and creates a dictionary of lineages based on indexes."""

def loadRankedLineage(indexes, rootNode):
    taxonomy = {}

    with open(ranked_lineage_filepath) as f:
        ranked_lineages = f.readlines()

    for lineage in ranked_lineages:
        lineage = extract_data(lineage)
        tax_index = lineage[0] # Put the taxonomy ID in a separate list.
        tax_index = "".join(tax_index.split()) # Clear messy spaces
        data = lineage[1:10] # take the actual ranked lineage into a separate list
        taxonomy[tax_index] = data # assign ranked lineage to its respective tax ID in order to a quickly accessible dictionary.

    for i in indexes:
        lineage = taxonomy[i]
        prepare_for_ranked_tree(lineage, rootNode)

"""Creates and linkes the nodes for the ETE Toolkit tree. Also applies some basic styling. """

def add_Edges(node, tree, rootNode):

    def calculate_thickness(child, style):

        thickness = (child.total_children / rootNode.total_children) * 60

        if child.rank == "Superkingdom" or child.rank == "Kingdom": # just to offset MASSIVE thick lines since those two
           thickness = 3                                            # will typically contain >90% of all results

        if thickness < 3: # just to offset later nodes having invisible borders
            thickness = 3

        style["vt_line_width"] = thickness
        style["hz_line_width"] = thickness
        return style #unecessary since it passes by reference, but it might be useful in future

    def highlight_tax(child, style):
        ranks = {'Name':'FireBrick',
                 'Species':'Crimson',
                 'Genus':'Chocolate',
                 'Family':'Gold',
                 'Order':'LawnGreen',
                 'Class':'LimeGreen',
                 'Phylum':'Turquoise',
                 'Kingdom':'SteelBlue',
                 'Superkingdom':'Plum'}
        style["bgcolor"] = ranks[child.rank]
        return style #unecessary since it passes by reference, but it might be useful in future


    for child in node.child_nodes:
        if child.rank is not ("Name"): # Restricting it to above Species means it won't crash rendering a tree with thousands of results

            newNode = Tree(name=child.node_name)

            # Set borders
            name = TextFace(child.node_name)
            name.border.width = 1
            name.border.color = '#ffffff'

            # Add name, total children, and rank.
            newNode.add_face(name,column=0,position="branch-right")
            children = "[" + str(child.total_children) + "]"
            newNode.add_face(TextFace(children),column=0, position="branch-bottom")
            newNode.add_face(TextFace(child.rank),column=0, position="branch-top")

            # Set colour based on rank and thickness based on total number of node's children.
            style = NodeStyle()
            highlight_tax(child, style)
            calculate_thickness(child, style)
            newNode.set_style(style)

            # Keep iterating
            tree.add_child(newNode)
            add_Edges(child, newNode, rootNode)

"""Creates the ETE tree needed for any sort of ETE Toolkit based tree visualisations. """

def create_ETE_tree(rootNode):
    rootETE = Tree()
    add_Edges(rootNode, rootETE, rootNode)

    return rootETE

def get_ASCII_tree(ETETree):
    return ETETree.get_ascii()

"""Prints the tree in terminal. """

def show_ASCII_tree(ETEtree):
        print(ETEtree.get_ascii())

"""Shows the tree in a GUI environment. """

def show_GUI_TREE(ETETree):
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_border = True

    ETETree.show(tree_style=ts)

"""Parses in the """

def get_most_least(argument, rootNode, descending):

    split_list = argument.split()

    while len(split_list) > 0:
        rank = split_list.pop(0)
        no_of_results = split_list.pop(0)
        rootNode.calculateCommonnessByRank(rank,no_of_results,rootNode,descending)

""""The main method of this file. Accepts a list of taxonomy indexes and the core program arguments. """

def main(indexes, args):


    rootNode = TaxNode("ROOT")
    loadRankedLineage(indexes, rootNode)
    rootNode.countChildren()

    ETETree = None

    if args.most:
        get_most_least(args.most, rootNode, True)

    if args.least:
        get_most_least(args.least, rootNode, False)

    # Show ASCII tree
    if args.ascii:
        if  ETETree is None:
            ETETree = create_ETE_tree(rootNode)
            show_ASCII_tree(ETETree)
        else:
            show_ASCII_tree(ETETree)

    # Save ASCII tree to file
    if args.save_ascii:
        if ETETree is None:
            ETETree = create_ETE_tree(rootNode)

        if args.out == args.save_ascii:
            file = open(args.save_ascii, "a")
        else:
            file = open(args.save_ascii, 'w')
        file.write(get_ASCII_tree(ETETree))
        file.close()

    # Show GUI tree
    if args.visualise:
        if ETETree is None:
            ETETree = create_ETE_tree(rootNode)
            show_GUI_TREE(ETETree)
        else:
            show_GUI_TREE(ETETree)
