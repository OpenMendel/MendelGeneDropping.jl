"""
This module simulates gene flow by gene dropping.
"""
module MendelGeneDropping
#
# Required OpenMendel packages and modules.
#
using MendelBase
# using DataStructures                  # Now in MendelBase.
# using GeneralUtilities                # Now in MendelBase.
#
# Required external modules.
#
using DataFrames                        # From package DataFrames.

export GeneDropping

"""
This is the wrapper function for the Gene Dropping analysis option.
"""
function GeneDropping(control_file = ""; args...)

  const GENE_DROPPING_VERSION :: VersionNumber = v"0.1.0"
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println("  Gene Dropping analysis option")
  println("        version ", GENE_DROPPING_VERSION)
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
    # Allowed values are: "Unordered", "Ordered", "Sourced", or "Population".
  keyword["gene_drop_output"] = "Unordered"
  keyword["interleaved"] = true
  keyword["keep_founder_genotypes"] = false
  keyword["missing_data_pattern"] = "ExistingData" # Not yet implemented.
  keyword["missing_rate"] = 0.0
  keyword["repetitions"] = 1
  #
  # Process the run-time user-specified keywords that will control the analysis.
  # This will also initialize the random number generator.
  #
  process_keywords!(keyword, control_file, args)
  #
  # Check that the correct analysis option was specified.
  #
  lc_analysis_option = lowercase(keyword["analysis_option"])
  if (lc_analysis_option != "" &&
      lc_analysis_option != "genedropping")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
  end
  keyword["analysis_option"] = "GeneDropping"
  #
  # Read the genetic data from the external files named in the keywords.
  #
  (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Execute the specified analysis.
  #
  println(" \nAnalyzing the data.\n")
  execution_error = false
  new_pedigree_frame = genedropping_option(pedigree, person,
    locus, locus_frame, phenotype_frame, pedigree_frame, keyword)
  show(new_pedigree_frame)
  if execution_error
    println(" \n \nERROR: Mendel terminated prematurely!\n")
  else
    println(" \n \nMendel's analysis is finished.\n")
  end
  #
  # Finish up by closing, and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing
end # function GeneDropping

"""
This function conducts gene dropping at the model loci.
The results are placed in a dataframe.
The process is controlled by the keywords gene_drop_output
(with values Unordered, Ordered, Sourced, and Population);
missing_rate; new_pedigree_file; repetitions; seed; and
and the boolean keyword 'interleaved'.
"""
function genedropping_option(pedigree::Pedigree, person::Person,
  locus::Locus, locus_frame::DataFrame, 
  phenotype_frame::DataFrame, pedigree_frame::DataFrame,
  keyword::Dict{AbstractString, Any})

  (rows, columns) = size(pedigree_frame)
  new_pedigree_frame = deepcopy(pedigree_frame[1, 1:columns])
  gene_drop_output = keyword["gene_drop_output"]
  #
  # Determine whether to use the ordered or unordered allele separator in output.
  # In either case use the first character in the defined strings.
  #
  if gene_drop_output == "Unordered"
    separator_char = keyword["allele_separator"][1]
  elseif gene_drop_output == "Ordered"
    separator_char = keyword["ordered_allele_separator"][1]
  elseif gene_drop_output == "Sourced"
    separator_char = keyword["ordered_allele_separator"][1]
  elseif gene_drop_output == "Population"
    separator_char = keyword["ordered_allele_separator"][1]
  else
    throw(ArgumentError(
      "Illegal value for gene_drop_output. OpenMendel terminated!\n \n"))
  end
  separator = ascii(string(separator_char))
  #
  # Loop over all pedigree/repeat combinations.
  #
  if keyword["interleaved"]
    for repetition = 1:keyword["repetitions"]
      for ped = 1:pedigree.pedigrees
        #
        # Copy the relevant part of the pedigree data frame.
        #
        start = pedigree.start[ped]
        finish = pedigree.twin_finish[ped]
        extent = finish - start + 1
        frame = deepcopy(pedigree_frame[start:finish, 1:columns])
        #
        # Perform gene dropping.
        #
        (sampled_genotype, source) =
          simulate_genotypes(pedigree, person, locus, keyword, ped)
        #
        # Name the pedigree replicate.
        #
        frame[1:extent, :Pedigree] = pedigree.name[ped] * string(repetition)
        #
        # Output can be sampled alleles or sources.
        #
        if gene_drop_output == "Unordered" || gene_drop_output == "Ordered"
          converted_genotype = convert_sampled_genotype(locus,
            sampled_genotype, separator, gene_drop_output)
        else
          converted_genotype = convert_sampled_genotype(locus,
            source, separator, gene_drop_output)
        end
        #
        # Load the simulated data.
        #
        for loc = 1:locus.loci
          l = locus.locus_field_in_pedigree_frame[loc]
          frame[1:extent, l] = converted_genotype[:, loc]
        end
        #
        # Sort the new pedigree so that it is the same order as the input file.
        # NB: Thus the new pedigree frame will not agree
        # with the order of the individuals in the data structures.
        # Append the simulated data to the new pedigree frame.
        #
        sort!(frame, cols = [:EntryOrder])
        append!(new_pedigree_frame, frame)
      end
    end
  else
    for ped = 1:pedigree.pedigrees
      #
      # Copy the relevant part of the pedigree data frame.
      #
      start = pedigree.start[ped]
      finish = pedigree.finish[ped]
      extent = finish - start + 1
      #
      # Perform gene dropping.
      #
      for repetition = 1:keyword["repetitions"]
        frame = deepcopy(pedigree_frame[start:finish, 1:columns])
        (sampled_genotype, source) =
          simulate_genotypes(pedigree, person, locus, keyword, ped)
        #
        # Name the pedigree replicate.
        #
        frame[1:extent, :Pedigree] = pedigree.name[ped] * string(repetition)
        #
        # Output can be sampled alleles or sources.
        #
        if gene_drop_output == "Unordered" || gene_drop_output == "Ordered"
          converted_genotype = convert_sampled_genotype(locus,
            sampled_genotype, separator, gene_drop_output)
        else
          converted_genotype = convert_sampled_genotype(locus,
            source, separator, gene_drop_output)
        end
        #
        # Load the simulated data.
        #
        for l = 1:locus.model_loci
          l = locus.locus_field_in_pedigree_frame[loc]
          frame[1:extent, l] = converted_genotype[:, loc]
        end
        #
        # Sort the new pedigree so that it is the same order as the input file.
        # NB: Thus the new pedigree frame will not agree
        # with the order of the individuals in the data structures.
        # Append the simulated data to the new pedigree frame.
        #
        sort!(frame, cols = [:EntryOrder])
        append!(new_pedigree_frame, frame)
      end
    end
  end
  #
  # If an output file was named for the simulated pedigrees,
  # then output the data, without the EntryOrder or Inverse_Perm fields.
  #
  deleterows!(new_pedigree_frame, 1)
  new_pedigree_file = keyword["new_pedigree_file"]
  if new_pedigree_file != ""
    names_list = names(new_pedigree_frame)
    deleteat!(names_list, findin(names_list, [:EntryOrder]))
    deleteat!(names_list, findin(names_list, [:Inverse_Perm]))
    writetable(new_pedigree_file, new_pedigree_frame[:, names_list])
  end
  return new_pedigree_frame
end # function gene_dropping_option

"""
This function carries out gene dropping for a single pedigree
over the model loci. These are the loci common to the Pedigree
and Locus frames.
"""
function simulate_genotypes(pedigree::Pedigree, person::Person,
  locus::Locus, keyword::Dict{AbstractString, Any}, ped::Int)

  theta = locus.theta
  #
  # Assign a source label to each founder gene.
  #
  source =
    founder_source(pedigree, person, locus, ped, keyword["gene_drop_output"])
  #
  # Compute pedigree attributes and offsets.
  #
  populations = person.populations
  founders = pedigree.founders[ped]
  ped_size = pedigree.individuals[ped]
  descendants = ped_size - founders
  model_loci = locus.model_loci
  q = pedigree.start[ped] - 1 # offset for founders
  r = q + founders # offset for descendants
  #
  # Allocate arrays.
  #
  sampled_genotype = zeros(Int, 2, ped_size, model_loci)
  gene = ones(Int, 2, ped_size)
  meiosis = ones(Int, 2, descendants)
  old_meiosis = ones(Int, 2, descendants)
  #
  # Loop over all model loci.
  #
  for l = 1:model_loci
    loc = locus.model_locus[l]
    xlinked = locus.xlinked[loc]
    #
    # Choose a compatible ordered genotype for each founder.
    # Each allele frequency is a convex combination of the
    # allele frequencies in the component populations.
    # The frequency variable is a column vector.
    #
    for i = 1:founders
      frequency =
        vec(transpose(transpose(person.admixture[i + q, :]) *
	              locus.frequency[loc]))
      male = person.male[i + q]
      #
      # Generate a founder genotype at the current locus.
      #
      if keyword["keep_founder_genotypes"]
        genotype_set = person.genotype[i + q, l]
        gene[:, i] = choose_genotype(frequency, genotype_set, xlinked, male)
      else
        gene[:, i] = random_genotype(frequency, xlinked, male)
      end
    end
    #
    # Propagate these choices to the descendants.
    #
    for i = 1:descendants
      for j = 1:2
        if j == 1
          parent = person.mother[i + r] - q
        else
          parent = person.father[i + r] - q
        end
        #
        # Determine the parental source of each gamete and store
        # a value of 1 or 2 in the meiosis indicator.
        #
        uniform = rand()
        if l == 1
          meiosis[j, i] = round(Integer, uniform) + 1
        else
          #
          # Take into account recombination after the first locus.
          #
          if uniform < theta[j, l - 1]
            meiosis[j, i] = 3 - old_meiosis[j, i]
          else
            meiosis[j, i] = old_meiosis[j, i]
          end
        end
        #
        # Record the allele passed and the source.
        #
        gene[j, i + founders] = gene[meiosis[j, i], parent]
        source[j, i + founders, l] = source[meiosis[j, i], parent, l]
      end
      if xlinked && person.male[i + r]
        source[2, i + founders, l] = source[1, i + founders, l]
      end
    end
    #
    # Save the sampled genotypes.
    #
    for j = 1:ped_size
      for i = 1:2
        sampled_genotype[i, j, l] = gene[i, j]
      end
    end
    #
    # Save the meiosis indicators.
    #
    if l < locus.model_loci
      old_meiosis = deepcopy(meiosis)
    end
  end
  #
  # Extend the simulation to co-twins.
  #
  for i = pedigree.finish[ped] + 1:pedigree.twin_finish[ped]
    m = person.primary_twin[i]
    sampled_genotype[:, i - q, :] = sampled_genotype[:, m - q, :]
    source[:, i - q, :] = source[:, m - q, :]
  end
  #
  # Delete genotypes that should be missing.
  #
  for l = 1:locus.model_loci
    for i = 1:ped_size
      if rand() < keyword["missing_rate"]
        sampled_genotype[:, i, l] = [0 0]
      end
    end
  end
  return (sampled_genotype, source)
end # function simulate_genotypes

"""
This function computes source labels for founder genes in a single pedigree.
Gene sources can be labeled by founders or ancestral populations.
"""
function founder_source(pedigree::Pedigree, person::Person,
  locus::Locus, ped::Int, gene_drop_output::AbstractString)

  q = pedigree.start[ped] - 1
  founders = pedigree.founders[ped]
  ped_size = pedigree.individuals[ped]
  model_loci = locus.model_loci
  source = zeros(Int, 2, ped_size, model_loci)
  #
  # Assign a source number to each founder gene.
  #
  for i = 1:founders
    if gene_drop_output == "Population"
      for l = 1:model_loci
        source[1, i, l] = random_category(vec(person.admixture[i + q, :]))
        source[2, i, l] = random_category(vec(person.admixture[i + q, :]))
      end
    else
      source[1, i, :] = 2i - 1
      source[2, i, :] = 2i
    end
  end
  #
  # At an x-linked locus, male founders have only one source.
  #
  for l = 1:model_loci
    loc = locus.model_locus[l]
    if locus.xlinked[loc]
      for i = 1:founders
        if person.male[i + q]
          source[2, i, l] = source[1, i, l]
        end
      end
    end
  end
  return source
end # function founder_source

"""
This function randomly chooses an ordered genotype from
a set of ordered genotypes. If the set is empty,
it returns a completely random genotype.
"""
function choose_genotype(frequency::Vector{Float64},
  genotype_set::Set{Tuple{Int, Int}}, xlinked::Bool, male::Bool)

  if length(genotype_set) == 0
    return random_genotype(frequency, xlinked, male)
  end
  #
  # Normalize the probabilities of the possible genotypes.
  #
  prob = zeros(length(genotype_set))
  i = 0
  for (gm, gp) in genotype_set
    i = i + 1
    if xlinked && male
      if gm == gp
        prob[i] = frequency[gm]
      end
    else
      prob[i] = frequency[gm] * frequency[gp]
    end
  end
  #
  # Check whether any compatible genotypes exist.
  #
  total = sum(prob)
  if total <= 0.0
    return random_genotype(frequency, xlinked, male)
  end
  #
  # Normalize the probabilities and return a sampled genotype.
  #
  prob = prob / total
  j = random_category(prob)
  i = 0
  for (gm, gp) in genotype_set
    i = i + 1
    if i == j
      return [gm gp]
    end
  end
end # function choose_genotype

"""
This function samples a random ordered genotype from
the universe of possible ordered genotypes.
"""
function random_genotype(frequency::Vector{Float64}, xlinked::Bool, male::Bool)

  if xlinked && male
    i = random_category(frequency)
    return [i i]
  else
    i = random_category(frequency)
    j = random_category(frequency)
    return [i j]
  end
end # function random_genotype

"""
This function converts sampled numerical genotypes into ordinary genotypes.
"""
function convert_sampled_genotype(locus::Locus, sampled_genotype::Array{Int, 3},
  separator::AbstractString, gene_drop_output::AbstractString)

  (m, n) = (size(sampled_genotype, 2), size(sampled_genotype, 3))
  converted_genotype = Array{AbstractString}(m, n)
  if gene_drop_output == "Unordered" || gene_drop_output == "Ordered"
    #
    # Convert allele numbers to allele names and concatenate.
    #
    for l = 1:locus.model_loci
      loc = locus.model_locus[l]
      for j = 1:m
        a = locus.allele_name[loc][sampled_genotype[1, j, l]]
        b = locus.allele_name[loc][sampled_genotype[2, j, l]]
        converted_genotype[j, l] = a * separator * b
      end
    end
  else
    #
    # Convert alleles numbers to integer strings and concatenate.
    #
    for l = 1:locus.model_loci
      loc = locus.model_locus[l]
      for j = 1:m
        a = string(sampled_genotype[1, j, l])
        b = string(sampled_genotype[2, j, l])
        converted_genotype[j, l] = a * separator * b
      end
    end
  end
  return converted_genotype
end # function convert_sampled_genotype

end # module MendelGeneDropping

