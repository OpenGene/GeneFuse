using OpenGene
using OpenGene.Reference
using ArgParse

function main(args)

    s = ArgParseSettings(description = "generate a fusion CSV file for given gene list")

    @add_arg_table s begin
        "--genes",  "-g"
            help = "the file giving gene list (separated by space, tab or line break)"
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

    passed_num = 0
    for genestr in allgenes
        genelist = split(genestr, '-')
        for genename in genelist
            if genename in processed_genes
                #println("----skipped ", genename, ", already processed")
                continue
            end
            gene_attrs = split(genename, '_')
            gene = gene_attrs[1]
            transcript = ""
            if length(gene_attrs)>1
                transcript = gene_attrs[2]
            end
            gene_matches = gencode_genes(index, gene)
            if length(gene_matches)==0
                #println("----skipped ", genename, ", gene name not found in gencode")
                continue
            end
            println(genename)
            g = gene_matches[1]

            # get the correct transcript if it is specified
            target_tr = nothing
            if transcript!=""
                for tr in g.transcripts
                    if startswith(tr.id, transcript)
                        target_tr = tr
                        break
                    end
                end
                if target_tr == nothing
                    #println("----skipped ", genename, ", transcript not found")
                end
            end

            # get the 001 transcript
            if target_tr == nothing
                for tr in g.transcripts
                    tname = tr.attributes["transcript_name"]
                    if endswith(tname, "001\"")
                        target_tr = tr
                        break
                    end
                end
            end

            # get the longest transcript if it is not specified
            if target_tr == nothing
                longest = 0
                for tr in g.transcripts
                    if length(tr.exons)>longest
                        longest = length(tr.exons)
                        target_tr = tr
                    end
                end
            end

            write(output, ">", gene, "_", target_tr.id, ",", g.chr, ":", string(g.start_pos), "-", string(g.end_pos), "\n")

            exons = target_tr.exons
            id = 1
            for e in exons
                write(output, string(id), ",", string(e.start_pos), ",", string(e.end_pos), "\n")
                id += 1
            end
            passed_num += 1
            write(output, "\n")
            push!(processed_genes, genename)
        end
    end
    println("successfully processed genes: ", passed_num)
    flush(output)
end

# REPL is not supported

main(ARGS)
