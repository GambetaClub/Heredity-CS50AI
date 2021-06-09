import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue
        
    # Uncomment the section below to stop testing 
        
        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def get_inheritance_prob(parent_genes, is_inherited):
    """
    Accepts the amount of genes the parent has and the boolean 
    descrbing whether to get the probability of inherit a gene or not. 
    it returns the probability of getting the gene (or not) from his parent.
    """
    if is_inherited:
        if parent_genes == 0:
            return PROBS["mutation"]
        elif parent_genes == 1:
            return 0.5
        else:
            return 1 - PROBS["mutation"]
    else:
        if parent_genes == 0:
            return 1 - PROBS["mutation"]
        elif parent_genes == 1:
            return 0.5
        else:
            return PROBS["mutation"]
    

def get_gene_quantity(person, one_gene, two_genes):
    if person in two_genes:
        return 2
    elif person in one_gene:
        return 1
    else:
        return 0

def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    
    probability = 1
    record = {}
    
    # For every person in people
    for person in people:

        calculation = 1
        # First, we see what we have to compute
        # See the probability of person of having the specified amount of genes
        gene_number = get_gene_quantity(person, one_gene, two_genes)
        
        # see the probability if person having the trait or not
        has_trait = person in have_trait
        

        # Easier accessors for father and mother 
        father = people[person]["father"]
        mother = people[person]["mother"]
        
        # Chech whether the person has father and mother
        if mother is None and father is None:
            # If it doesn't, then compute the inconditional probability
            calculation *= PROBS["gene"][gene_number] * PROBS["trait"][gene_number][has_trait]
        else:
            father_genes = get_gene_quantity(father, one_gene, two_genes)
            mother_genes = get_gene_quantity(mother, one_gene, two_genes)
            # If it does, then compute the probability depending on their parents
            calculation = 0
            #  Probability of person for getting two genes
            if gene_number == 2:
                # Only if they inherit both from their father and mother
                calculation += get_inheritance_prob(father_genes, True) * get_inheritance_prob(father_genes, True)
            #  Probabilities of person for getting only one gene
            elif gene_number == 1:
                # If they inherit from their father and not their mother or viceversa
                calculation += get_inheritance_prob(father_genes, True) * get_inheritance_prob(father_genes, False)
                calculation += get_inheritance_prob(father_genes, False) * get_inheritance_prob(father_genes, True)
            #  Probabilities of person for getting zero genes
            else:
                
                # Only if they don't inherit from neither from their nor their mother 
                 calculation += get_inheritance_prob(father_genes, False) * get_inheritance_prob(father_genes, False)
            
        record[person] = {
        "Name": person,
        "Number of Genes": gene_number,
        "Trait": has_trait,
        "Probability": calculation,
        }
        
        probability *= calculation
    
    # for key, value in record.items():
    #     print(key, ' : ', value)
    
    return probability

def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for person in probabilities:
        gene_number = get_gene_quantity(person, one_gene, two_genes)
        has_trait = person in have_trait
        probabilities[person]["gene"][gene_number] += p
        probabilities[person]["trait"][has_trait] += p


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in probabilities:
        trait_a = sum(probabilities[person]["trait"].values())
        gene_a = sum(probabilities[person]["gene"].values())
        
        for gene in probabilities[person]["gene"]:
            probabilities[person]["gene"][gene] /= gene_a
            
        for trait in probabilities[person]["trait"]:
            probabilities[person]["trait"][trait] /= trait_a

    
            
if __name__ == "__main__":
    main()
