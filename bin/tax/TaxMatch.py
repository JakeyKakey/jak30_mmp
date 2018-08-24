import time
import pickle
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

"""Loads the provided taxonomy index/names file and splits it into two separate lists """

names_and_indexes_filepath = '../bin/database/taxonomy/names-trimmed.txt'

def load_taxonomy():
    tax_species = []
    tax_index = []

    start = time.time()

    with open(names_and_indexes_filepath) as f:
        index_matches = f.readlines()

        for species in index_matches:
            species = species.split("|")
            tax_species.append(species[1].lower())
            tax_index.append(species[0])

    end = time.time()
    print("Taxonomy list loaded in ", end - start)

    return tax_index, tax_species

"""Creates a bag-of-words TF-IDF matrix/model from the list of taxonomy species names """

def compute_tax_matrix(tax_species):

    print("Computing taxonomy matrix.")
    start = time.time()

    tfidf_vectorizer = TfidfVectorizer()
    matrix = tfidf_vectorizer.fit_transform(tax_species)

    end = time.time()

    print("Taxonomy matrix computed in ", end - start)

    return matrix, tfidf_vectorizer

"""Cleans up the list of subject names to remove special characters and some of the common stop words, e.g: 
    From: gi|551368999|emb|BX284682.8| Zebrafish DNA sequence from clone DKEYP-74C5 in linkage group 22, complete sequence
    To:   gi551368999embBX284682 8 Zebrafish DKEYP-74C5 22
"""

def format_results(results):

    # Mildly redundant as the vectoriser probably has a way to insert stop words,
    # but frankly I don't want to dive into that documentation anymore. Possibly refactor later.

    stop_words = ("complete sequence", " complete ", " sequence ", " scaffold ", " chromosome ",
                  " DNA ", " group ", " linkage ", " assembly ", " strain ",
                  " PREDICTED ", " hypothetical ", " partial ", " protein ",
                  " from ", " transcript ", " variant ", " clone ", " in ", " genome ", " contig ", " mRNA ")

    cleaned=[]

    for x in results:
        #print(x) #DEBUG
        x = x + " "
        x = x.replace("_"," ")
        x = x.replace(","," ")
        x = x.replace("|", "") # Intentionally replaced with no space to reduce the effect it will have on bag of words.
        x = x.replace(".", " ")
        x = x.replace("=", " ")
        x = x.replace(";", " ")
        x = x.replace(":", " ")

        for word in stop_words:
            x = x.replace(word, " ")
        cleaned.append(x)

    return cleaned

"""Creates a bag-of-words matrix to compare against the taxonomy matrix. 
    Requires the original vectorizer used in compute_taxonomy_matrix in order to work.
"""

def compute_matches_matrix(results, vectorizer):

    print("Computing matches matrix.")
    start = time.time()
    matches_matrix = vectorizer.transform(results)
    end = time.time()
    print("Matches matrix computed in ", end - start)

    return matches_matrix

"""Combines the bag-of-words term indexes and their cosine vectors into a regular list,
    sorts it by top results, and throws away the rest.
"""

def analyse_match(match):
    linked = {}
    import operator
    for x, y in zip(match.data, match.indices):
        linked[y] = x
    sorted_linked = sorted(linked.items(), key=operator.itemgetter(1), reverse=True)
    if len(sorted_linked) > 5: # change here and below if want to maintain more than 5 top results
        sorted_linked = sorted_linked[0:5]
    return sorted_linked

""" Very crude accuracy check. Tests if the matched species string is contained within the string of original BLAST match name.
    Real accuracy should be a bit higher as this is quite prone to false-negatives.
"""

def check_accuracy(original, matched):
    a = original.lower()
    b = matched.lower()
    b = " ".join(b.split())

    if b not in a:
        print("WRONG?")
        return True
    return False


def match_results_to_taxonomy(BLAST_matches):

    index_list = []

    clean_matches = format_results(BLAST_matches) # Clean up the results for bag-of-words.
    tax_index, tax_species = load_taxonomy() # Load the species names and their indexes.
    taxonomy_matrix, vectorizer = compute_tax_matrix(tax_species) #Create bag-of-words taxonomy matrix.
    matches_matrix = compute_matches_matrix(clean_matches, vectorizer) #Create bag-of-words species matrix.

    matches_size = matches_matrix.shape[0] # The direct array manipulation here as well as in analyse_match below,
                                           # is kind of a bad practice, but these are fundamentally just 2D arrays,
                                           # while Numpy multi-dimensional matrices are pretty confusing to operate on.


    mismatched = 0

    for i in range(0, matches_size):
        result = cosine_similarity(matches_matrix[i], taxonomy_matrix, False)
        likely_matches = analyse_match(result) # Will return the top 5 highest matches.



    #For debugging and printing out the top 5 matches.
        print((i+1),"out of:", matches_size, " Here's a list of likely matches for", BLAST_matches[i])
        for x in likely_matches:
            match_id = x[0]
            percentage = x[1]
            print(tax_index[match_id], tax_species[match_id], percentage)



    # Otherwise grab the top one, for better or worse.

        top_match = likely_matches[0]
        match_id = top_match[0]
        percentage = top_match[1]
        #print("Matched:", i,"out of", matches_size,".")
        #print(matches[i], "matched to", tax_species[match_id], "with", (percentage*100), "accuracy.") #DEBUG
        index = tax_index[match_id]
        index = "".join(index.split())
        index_list.append(index)

        if check_accuracy(BLAST_matches[i], tax_species[match_id]):
            mismatched += 1

    accuracy = 100 - ((mismatched / matches_size)*100)

    print("Matched with", accuracy, "% accuracy.")

    return index_list
