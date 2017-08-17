using OpenGene
using OpenGene.Reference
using ArgParse

function main(args)

    s = ArgParseSettings(description = "gene a fusion CSV file for given gene list")

    @add_arg_table s begin
        "--genes",  "-g"
            help = "the file giving gene list (separated by comma, space or line break)"
            required = true
        "--ref",  "-r"
            help = "the reference genome to use (hg19/hg37/grch38)"
            required = true
        "--fusion",   "-f"
            help = "the fusion file to output"
            required = true
    end

    options = parse_args(s) # the result is a Dict{String,Any}
    detect(options["ref"], options["genes"], options["fusion"])
end

function detect(ref, genes, fusion)
    ref = lowercase(ref)
    if !(ref in ["hg19", "grch37", "grch38"])
        println("reference <-f/--ref> should be one of hg19/grch37/grch38")
        return
    end
    input = open(genes)
    output = open(fusion, "w")
    genestr = readstring(input)
    allgenes = split(genestr)

    println("loading gencode " * ref)
    index = gencode_load(ref)
    println("loading done")

    processed_genes = []

    for genestr in allgenes
        genelist = split(genestr, ',')
        for gene in genelist
            gene = uppercase(gene)
            gene_matches = gencode_genes(index, gene)
            if length(gene_matches)==0
                println("----skipped ", gene, ", gene name not found in gencode")
                continue
            end
            if gene in processed_genes
                println("----skipped ", gene, ", duplicated")
                continue
            end
            println(gene)
            g = gene_matches[1]
            write(output, ">", g.name, ",", g.chr, ":", string(g.start_pos), "-", string(g.end_pos), "\n")
            exons = g.transcripts[1].exons
            id = 1
            for e in exons
                write(output, string(id), ",", string(e.start_pos), ",", string(e.end_pos), "\n")
                id += 1
            end
            write(output, "\n")
            push!(processed_genes, gene)
        end
    end
    flush(output)
end

# REPL is not supported

main(ARGS)
