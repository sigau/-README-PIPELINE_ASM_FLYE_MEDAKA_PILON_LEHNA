from multiprocessing import cpu_count
full_core=cpu_count()

configfile: "config.yaml"


rule porechop_reads:
    input:
        config["reads_folder"]
    output:
        temp(expand("{out_dir}/reads/all_reads_porechoped.fastq.gz", out_dir=config["out_dir"]))
    threads: full_core
    log: expand("{out_dir}/logs/porechop_reads/porechop.log", out_dir=config["out_dir"])
    message: "removing the adapter from the reads usig porechop"
    conda:
        "env/porechop.yaml"
    shell:
        "( porechop --no_split -t {threads} -i {input} -o {output} ) 2> {log}"


rule filter_reads:
    input:
        rules.porechop_reads.output
    output:
        temp(expand("{out_dir}/reads/filtered_reads_porechoped.fastq.gz", out_dir=config["out_dir"]))
    log: expand("{out_dir}/logs/filtering_reads/nanofilt.log", out_dir=config["out_dir"])
    message: "filtering the reads under min_len base and of quality inferior to min_q"
    params:
        q=config["min_q_LR"],
        l=config["min_len_LR"]
    conda:
        "env/nanofilt.yaml"
    shell:
        "( gunzip -c {input} | NanoFilt  -l {params.l} -q {params.q} | gzip > {output} ) 2> {log}"


rule QC_reads:
    input:
        rules.filter_reads.output
    output:
        expand("{out_dir}/QC/nanoplot/NanoPlot-report.html", out_dir=config["out_dir"])
    threads: 4
    log: expand("{out_dir}/logs/filtering_reads/nanofilt.log", out_dir=config["out_dir"])
    message: "QC of the filtering using nanoplot"
    conda:
        "env/nanofilt.yaml"
    params:
        folder=expand("{out_dir}/QC/nanoplot/", out_dir=config["out_dir"])
    shell:
        "( NanoPlot -t {threads} -o {params.folder} --fastq {input} ) 2> {log}"


rule assembly_flye:
    input:
        rules.filter_reads.output
    output:
        expand("{out_dir}/assembly_flye_hq/{asm_name}_asm_flyeHQ.fasta", out_dir=config["out_dir"], asm_name=config["asm_name"])
    threads: full_core
    log: expand("{out_dir}/logs/assembly/flye.log", out_dir=config["out_dir"])
    message: "Assembling the reads with flye using HQ setting"
    conda:
        "env/flye.yaml"
    params:
        gen_size=config["genome_size"],
        out_dir=expand("{out_dir}/assembly_flye_hq", out_dir=config["out_dir"])
    shell:
        "( flye --threads {threads} --nano-hq {input} --out-dir {params.out_dir} --genome-size {params.gen_size} ) 2> {log} && mv {params.out_dir}/assembly.fasta {output}"


rule polishing_medaka:
    input:
        asm=rules.assembly_flye.output,
        reads=rules.filter_reads.output
    output:
        expand("{out_dir}/polishing_medaka/{asm_name}_flyeHQ_medaka.fasta", out_dir=config["out_dir"], asm_name=config["asm_name"])
    threads: full_core
    message: "Polishing the assembly with medaka"
    params:
        model=config["medaka_model"],
        out_dir=expand("{out_dir}/polishing_medaka", out_dir=config["out_dir"])
    log: expand("{out_dir}/logs/polishing/mdeka.log", out_dir=config["out_dir"])
    conda:
        "env/medaka.yaml"
    shell:
        "(medaka_consensus -i {input.reads} -d {input.asm} -o {params.out_dir} -m {params.model} && cp {params.out_dir}/consensus.fasta {output} ) 2> {log}"

rule curing:
    input:
        R1=config["short_r1"],
        R2=config["short_r2"]
    output:
        r1=expand("{out_dir}/reads/R1_trim_paired.fastq.gz", out_dir=config["out_dir"]),
        r2=expand("{out_dir}/reads/R2_trim_paired.fastq.gz", out_dir=config["out_dir"]),
        r1u=expand("{out_dir}/reads/R1_trim_unpaired.fastq.gz", out_dir=config["out_dir"]),
        r2u=expand("{out_dir}/reads/R2_trim_unpaired.fastq.gz", out_dir=config["out_dir"])
    threads: full_core
    log: expand("{out_dir}/logs/Trimo/trimo.log", out_dir=config["out_dir"])
    message: "triming the short-reads with trimmomatic"
    params:
        min_len=config["min_len"],
        min_q=config["min_q"]
    shell:
        " ( /usr/local/jre1.8.0_202/bin/java -jar src/trimmomatic-0.39.jar PE -phred33 {input.R1} {input.R2} {output.r1} {output.r1u} "
        " {output.r2} {output.r2u} SLIDINGWINDOW:4:{params.min_q} MINLEN:{params.min_len} ) 2> {log}"


rule pilon1_bowtie:
    input:
        r1=rules.curing.output.r1,
        r2=rules.curing.output.r2,
        asm=rules.polishing_medaka.output
    output:
        sam=temp(expand("{out_dir}/bowtie/mapping_1.sam", out_dir=config["out_dir"])),
        bam=temp(expand("{out_dir}/bowtie/mapping_1.bam", out_dir=config["out_dir"])),
        bam_sort=temp(expand("{out_dir}/bowtie/mapping_1_sort.bam", out_dir=config["out_dir"]))
    threads: full_core
    log: expand("{out_dir}/logs/bowtie/bowtie1.log", out_dir=config["out_dir"])
    message: "mapping the reads against the assembly"
    params:
        index=expand("{out_dir}/bowtie/index_1", out_dir=config["out_dir"])
    conda:
        "env/bowtie2.yaml"
    shell:
        "( bowtie2-build {input.asm} {params.index} && "
        " bowtie2 --threads {threads} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.sam} && "
        " samtools view -bS {output.sam} > {output.bam} && samtools sort -@ {threads} -o {output.bam_sort}  {output.bam} &&"
        " samtools index {output.bam_sort} && bwa index {input.asm} ) 2> {log}"

rule pilon1_polish:
    input:
        asm=rules.polishing_medaka.output,
        bam=rules.pilon1_bowtie.output.bam_sort
    output:
        temp(expand("{out_dir}/pilon/round1/{asm_name}_pilon_round_1.fasta", out_dir=config["out_dir"], asm_name=config["asm_name"]))
    threads: full_core
    log: expand("{out_dir}/logs/pilon/pilon_1.log", out_dir=config["out_dir"])
    message: "polishing the assembly with pilon"
    params:
        output=expand("{out_dir}/pilon/round1/{asm_name}_pilon_round_1", out_dir=config["out_dir"], asm_name=config["asm_name"])
    shell:
        "( /usr/local/jre1.8.0_202/bin/java -Xmx64G -jar src/pilon-1.24.jar --threads {threads} --genome {input.asm} --frags {input.bam} --output {params.output} ) 2> {log}"

rule pilon2_bowtie:
    input:
        r1=rules.curing.output.r1,
        r2=rules.curing.output.r2,
        asm=rules.pilon1_polish.output
    output:
        sam=temp(expand("{out_dir}/bowtie/mapping_2.sam", out_dir=config["out_dir"])),
        bam=temp(expand("{out_dir}/bowtie/mapping_2.bam", out_dir=config["out_dir"])),
        bam_sort=temp(expand("{out_dir}/bowtie/mapping_2_sort.bam", out_dir=config["out_dir"]))
    threads: full_core
    log: expand("{out_dir}/logs/bowtie/bowtie2.log", out_dir=config["out_dir"])
    message: "mapping the reads against the assembly"
    params:
        index=expand("{out_dir}/bowtie/index_2", out_dir=config["out_dir"])
    conda:
        "env/bowtie2.yaml"
    shell:
        "( bowtie2-build {input.asm} {params.index} && "
        " bowtie2 --threads {threads} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.sam} && "
        " samtools view -bS {output.sam} > {output.bam} && samtools sort -@ {threads} -o {output.bam_sort}  {output.bam} &&"
        " samtools index {output.bam_sort} && bwa index {input.asm} ) 2> {log}"

rule pilon2_polish:
    input:
        asm=rules.pilon1_polish.output,
        bam=rules.pilon2_bowtie.output.bam_sort
    output:
        temp(expand("{out_dir}/pilon/round2/{asm_name}_pilon_round_2.fasta", out_dir=config["out_dir"], asm_name=config["asm_name"]))
    threads: full_core
    log: expand("{out_dir}/logs/pilon/pilon_2.log", out_dir=config["out_dir"])
    message: "polishing the assembly with pilon"
    params:
        output=expand("{out_dir}/pilon/round2/{asm_name}_pilon_round_2", out_dir=config["out_dir"], asm_name=config["asm_name"])
    shell:
        "( /usr/local/jre1.8.0_202/bin/java -Xmx64G -jar src/pilon-1.24.jar --threads {threads} --genome {input.asm} --frags {input.bam} --output {params.output} ) 2> {log}"

rule pilon3_bowtie:
    input:
        r1=rules.curing.output.r1,
        r2=rules.curing.output.r1,
        asm=rules.pilon2_polish.output
    output:
        sam=temp(expand("{out_dir}/bowtie/mapping_3.sam", out_dir=config["out_dir"])),
        bam=temp(expand("{out_dir}/bowtie/mapping_3.bam", out_dir=config["out_dir"])),
        bam_sort=temp(expand("{out_dir}/bowtie/mapping_3_sort.bam", out_dir=config["out_dir"]))
    threads: full_core
    log: expand("{out_dir}/logs/bowtie/bowtie3.log", out_dir=config["out_dir"])
    message: "mapping the reads against the assembly"
    params:
        index=expand("{out_dir}/bowtie/index_3", out_dir=config["out_dir"])
    conda:
        "env/bowtie2.yaml"
    shell:
        "( bowtie2-build {input.asm} {params.index} && "
        " bowtie2 --threads {threads} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.sam} && "
        " samtools view -bS {output.sam} > {output.bam} && samtools sort -@ {threads} -o {output.bam_sort}  {output.bam} &&"
        " samtools index {output.bam_sort} && bwa index {input.asm} ) 2> {log}"

rule pilon_polish:
    input:
        asm=rules.pilon2_polish.output,
        bam=rules.pilon3_bowtie.output.bam_sort
    output:
        expand("{out_dir}/pilon/round3/{asm_name}_medaka_pilon_round_3.fasta", out_dir=config["out_dir"], asm_name=config["asm_name"])
    threads: full_core
    log: expand("{out_dir}/logs/pilon/pilon_3.log", out_dir=config["out_dir"])
    message: "polishing the assembly with pilon"
    params:
        output=expand("{out_dir}/pilon/round3/{asm_name}_medaka_pilon_round_3", out_dir=config["out_dir"], asm_name=config["asm_name"])
    shell:
        "( /usr/local/jre1.8.0_202/bin/java -Xmx64G -jar src/pilon-1.24.jar --threads {threads} --genome {input.asm} --frags {input.bam} --output {params.output} ) 2> {log}"    

rule QC_assembly_quast:
    input:
        rules.assembly_flye.output
    output:
        expand("{out_dir}/QC/QUAST/DRAFT_ASSEMBLY/report.tsv", out_dir=config["out_dir"])
    threads: 4 
    log: expand("{out_dir}/logs/QC_QUAST/draft_assembly.log", out_dir=config["out_dir"])
    message: "Quality Control of the assembly using QUAST"
    params:
        out_dir=expand("{out_dir}/QC/QUAST/DRAFT_ASSEMBLY", out_dir=config["out_dir"])
    conda:
        "env/quast.yaml"
    shell:
        "(quast {input} -o {params.out_dir} -t {threads} --eukaryote --large) 2> {log}"


rule QC_assembly_busco:
    input:
        rules.assembly_flye.output
    output:
        expand("{out_dir}/QC/BUSCO/{asm_name}_DRAFT/logs/busco.log", out_dir=config["out_dir"],asm_name=config["asm_name"])
    threads: full_core
    message: "Quality Control of the assembly using BUSCO"
    params:
        out_dir=expand("{out_dir}/QC/BUSCO/", out_dir=config["out_dir"]),
        busco_name=expand("{asm_name}_DRAFT",asm_name=config["asm_name"]),
        db=config["busco_db"]
    log: expand("{out_dir}/logs/QC_BUSCO/draft_assembly.log", out_dir=config["out_dir"])
    conda:
        "env/busco.yaml"
    shell:
        "( busco -f -c {threads} -l {params.db} -m genome --out_path {params.out_dir} -i {input} -o  {params.busco_name}) 2> {log}"


rule QC_polishing_quast_SR:
    input:
        rules.pilon_polish.output
    output:
        expand("{out_dir}/QC/QUAST/POLISHING_ASSEMBLY_MEDAKA_PILON/report.tsv", out_dir=config["out_dir"])
    threads: 4 
    log: expand("{out_dir}/logs/QC_QUAST/polish_assembly.log", out_dir=config["out_dir"])
    message: "Quality Control of the polish assembly using QUAST"
    params:
        out_dir=expand("{out_dir}/QC/QUAST/POLISHING_ASSEMBLY_MEDAKA_PILON", out_dir=config["out_dir"])
    conda:
        "env/quast.yaml"
    shell:
        "(quast {input} -o {params.out_dir} -t {threads} --eukaryote --large) 2> {log}"


rule QC_polishing_busco_SR:
    input:
        rules.pilon_polish.output
    output:
        expand("{out_dir}/QC/BUSCO/{asm_name}_medaka_pilon_polish/logs/busco.log", out_dir=config["out_dir"],asm_name=config["asm_name"])
    threads: full_core
    message: "Quality Control of the polish assembly using BUSCO"
    params:
        out_dir=expand("{out_dir}/QC/BUSCO/", out_dir=config["out_dir"]),
        busco_name=expand("{asm_name}_medaka_pilon_polish",asm_name=config["asm_name"]),
        db=config["busco_db"]
    log: expand("{out_dir}/logs/QC_BUSCO/polish_assembly.log", out_dir=config["out_dir"])
    conda:
        "env/busco.yaml"
    shell:
        "( busco -f -c {threads} -l {params.db} -m genome --out_path {params.out_dir} -i {input} -o {params.busco_name} ) 2> {log}"


rule QC_polishing_quast_LR:
    input:
        rules.polishing_medaka.output
    output:
        expand("{out_dir}/QC/QUAST/POLISHING_ASSEMBLY_MEDAKA/report.tsv", out_dir=config["out_dir"])
    threads: 4 
    log: expand("{out_dir}/logs/QC_QUAST/polish_assembly.log", out_dir=config["out_dir"])
    message: "Quality Control of the polish assembly using QUAST"
    params:
        out_dir=expand("{out_dir}/QC/QUAST/POLISHING_ASSEMBLY_MEDAKA", out_dir=config["out_dir"])
    conda:
        "env/quast.yaml"
    shell:
        "(quast {input} -o {params.out_dir} -t {threads} --eukaryote --large) 2> {log}"


rule QC_polishing_busco_LR:
    input:
        rules.polishing_medaka.output
    output:
        expand("{out_dir}/QC/BUSCO/{asm_name}_medaka_polish/logs/busco.log", out_dir=config["out_dir"],asm_name=config["asm_name"])
    threads: full_core
    message: "Quality Control of the polish assembly using BUSCO"
    params:
        out_dir=expand("{out_dir}/QC/BUSCO/", out_dir=config["out_dir"]),
        busco_name=expand("{asm_name}_medaka_polish",asm_name=config["asm_name"]),
        db=config["busco_db"]
    log: expand("{out_dir}/logs/QC_BUSCO/polish_assembly.log", out_dir=config["out_dir"])
    conda:
        "env/busco.yaml"
    shell:
        "( busco -f -c {threads} -l {params.db} -m genome --out_path {params.out_dir} -i {input} -o {params.busco_name} ) 2> {log}"


rule asm_polish_SR_LR:
    input:
        qp=rules.QC_polishing_quast_SR.output,
        bp=rules.QC_polishing_busco_SR.output,
        qa=rules.QC_assembly_quast.output,
        ba=rules.QC_assembly_busco.output,
        qp2=rules.QC_polishing_quast_LR.output,
        bp2=rules.QC_polishing_busco_LR.output,
        qr=rules.QC_reads.output

rule asm_polish_LR:
    input:
        qp=rules.QC_polishing_quast_LR.output,
        bp=rules.QC_polishing_busco_LR.output,
        qa=rules.QC_assembly_quast.output,
        ba=rules.QC_assembly_busco.output,
        qr=rules.QC_reads.output

rule assembly:
    input:
        qa=rules.QC_assembly_quast.output,
        ba=rules.QC_assembly_busco.output,
        qr=rules.QC_reads.output
