import os
from pathlib import Path

configfile: "configs/variant_calling.yaml"
configfile: "configs/variant_calling2.yaml"

ref = config["ref"]
outdir = config["outdir"]
example_name = config["example_name"]
vc_env = config["vc_env"]
vc_env2 = config["vc_env2"]
vc_env3 = config["vc_env3"]
vc_tools=config["vc_tools"]

def get_sim_config():
    configs = []
    no_indel = True
    for paras in config["simuG_paras"]:
        paras = paras.split(" ")
        assert len(paras) == 4
        if int(paras[3]) == 0:
            no_indel = True
            configs.append({
                "seed": paras[0],
                "ratio": paras[1],
                "snp_count": paras[2],             
            })

        else:
            no_indel = False
            configs.append({
                "seed": paras[0],
                "ratio": paras[1],
                "snp_count": paras[2],
                "indel_count": paras[3]               
            })
    return configs, no_indel

def get_simulation_files(example_name):
    configs, no_indel = get_sim_config()
    
    files = []
    fa_files = []
    for cfg in configs:
        files.append(
            f"{outdir}/simulated/{example_name}_seed{cfg['seed']}.refseq2simseq.SNP"
        )
        if not no_indel:
            files.append(
                f"{outdir}/simulated/{example_name}_seed{cfg['seed']}.refseq2simseq.INDEL"
            )
        fa_files.append(f"{outdir}/simulated/{example_name}_seed{cfg['seed']}.simseq.genome.fa")
    
    return files, fa_files

vcfs_name, fa_files = get_simulation_files(example_name)
snp_vcfs_name = [name.split("/")[-1].split(".")[0] for name in vcfs_name if "SNP" in name]

rule all:
    input:
        f"{outdir}/reference/ref.fa.fai",
        f"{outdir}/reference/ref.dict",
        f"{outdir}/reference/ref.sdf",
        expand("{vcf_name}.vcf", vcf_name=vcfs_name),
        expand("{vcf_name}.seq.vcf", vcf_name=vcfs_name),
        f"{outdir}/sim_reads/config.ini",
        f"{outdir}/sim_reads/camisim.log",
        f"{outdir}/simulated/{example_name}_variant_all.vcf.gz",
        f"{outdir}/simulated/{example_name}_variant_all.alt2.vcf.gz",
        f"{outdir}/align/{example_name}_align.bam",
        expand(outdir + "/tools_result/{vc_tool}/" + example_name + ".vcf.gz", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/eval.txt", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/rtg_eval.txt", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/rtg_eval_alt2.txt", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/rtg_eval_alt3.txt", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/hap_py_eval.txt", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/rtg_eval_{vc_tool}_{vcf_name_single}.txt", vc_tool=vc_tools, vcf_name_single=snp_vcfs_name),
        expand(outdir + "/tools_result/{vc_tool}/my_eval_{vc_tool}_{vcf_name_single}.txt", vc_tool=vc_tools, vcf_name_single=snp_vcfs_name),
        f"{outdir}/all_eval.md"
        
rule copy_ref:
    input:
        ref=ref,
    output:
        ref_fa=f"{outdir}/reference/ref.fa",
    shell:
        "mkdir -p {outdir}/reference && cp {input.ref} {output}"


rule process_ref:
    input:
        ref=f"{outdir}/reference/ref.fa",
    output:
        directory(f"{outdir}/reference/ref.sdf"),
        f"{outdir}/reference/ref.fa.fai",
        f"{outdir}/reference/ref.dict",
    conda:
        vc_env
    shell:
        """
        mkdir -p {outdir}/reference && cd {outdir}/reference
        samtools faidx ref.fa
        rtg format -o ref.sdf ref.fa
        gatk CreateSequenceDictionary -R ref.fa
        """





rule simulate_single:
    input:
        ref_fa=rules.copy_ref.output.ref_fa,
    output:
        expand("{vcf_name}.vcf", vcf_name=vcfs_name)
    params:
        simuG=config["simuG"]
    conda:
        vc_env
    run:
        import subprocess
    
        configs, no_indel = get_sim_config()
        
        for cfg in configs:
            seed = cfg['seed']
            raw_vcf = f"{outdir}/simulated/{example_name}_seed{seed}.refseq2simseq.SNP.vcf"
            seq_vcf = f"{outdir}/simulated/{example_name}_seed{seed}.refseq2simseq.SNP.seq.vcf"
            if not Path(raw_vcf).exists():
                if no_indel:
                    cmd = [
                        "perl", params.simuG,
                        "-refseq", input.ref_fa,
                        "-snp_count", str(cfg["snp_count"]),
                        "-titv_ratio", str(cfg["ratio"]),
                        "-seed", str(cfg["seed"]),
                        "-prefix", f"{outdir}/simulated/{example_name}_seed{cfg['seed']}"
                    ]
                else:
                    cmd = [
                        "perl", params.simuG,
                        "-refseq", input.ref_fa,
                        "-snp_count", str(cfg["snp_count"]),
                        "-indel_count", str(cfg["indel_count"]),
                        "-titv_ratio", str(cfg["ratio"]),
                        "-seed", str(cfg["seed"]),
                        "-prefix", f"{outdir}/simulated/{example_name}_seed{cfg['seed']}"
                    ]
                subprocess.run(cmd, check=True)

rule update_vcf:
    input: 
        vcf = "{vcf_name}.vcf",
        ref = rules.copy_ref.output.ref_fa
    output:
        seq_vcf = "{vcf_name}.seq.vcf"
    conda:
        vc_env
    shell:
        """
        gatk UpdateVCFSequenceDictionary -V {input.vcf} -O {output.seq_vcf} -R {input.ref}
        """

rule merge_all_vcfs:
    input:
        vcfs=expand("{vcf_name}.seq.vcf", vcf_name=vcfs_name)
    output:
        merged=f"{outdir}/simulated/{example_name}_merged.vcf.gz"
    conda:
        vc_env
    params:
        input_args = lambda wildcards, input: ' '.join([f'-I {f}' for f in input.vcfs])
    shell:
        """
        gatk MergeVcfs {params.input_args} -O {output.merged} 
        """

rule process_vcf:
    input:
        merged_vcf=rules.merge_all_vcfs.output.merged,
    output:
        final_vcf=f"{outdir}/simulated/{example_name}_variant_all.vcf.gz",
        final_tbi=f"{outdir}/simulated/{example_name}_variant_all.vcf.gz.tbi"
    conda:
        vc_env
    shell:
        """
        # set FORMAT and SAMPLE columns, need complete vcf format
        bcftools view {input.merged_vcf} | \
        awk 'BEGIN{{OFS="\t"}} /^#/ {{if ($0 ~ /^#CHROM/) print $0, "FORMAT", "SAMPLE"; else print $0; next}} {{print $0, "GT", "1/1"}}' | \
        bcftools annotate -h <(echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">') - | \
        bcftools sort -Ou | \
        bcftools norm -m +both -Oz -o {outdir}/simulated/{example_name}_variant_all.vcf.gz 
        tabix -p vcf {output.final_vcf}       
        """        

rule filter_alt_counts:
    input:
        vcf=rules.process_vcf.output.final_vcf,
        tbi=rules.process_vcf.output.final_tbi,
    output:
        alt2_vcf=f"{outdir}/simulated/{example_name}_variant_all.alt2.vcf.gz",
        alt2_tbi=f"{outdir}/simulated/{example_name}_variant_all.alt2.vcf.gz.tbi",
        alt3_vcf=f"{outdir}/simulated/{example_name}_variant_all.alt3.vcf.gz",
        alt3_tbi=f"{outdir}/simulated/{example_name}_variant_all.alt3.vcf.gz.tbi"
    conda:
        vc_env
    shell:
        """
        # filter 3 alts
        bcftools view -i 'N_ALT=3' {input.vcf} -o {output.alt3_vcf}
        tabix -p vcf {output.alt3_vcf}
        
        # filter 2 alts
        bcftools view -i 'N_ALT=2' {input.vcf} -o {output.alt2_vcf}
        tabix -p vcf {output.alt2_vcf}
        """

def create_camisim_config(sim_dir, fa_files, distribution, config_ini, size):
    sim_dir_path = Path(sim_dir)
    sim_dir_path.mkdir(parents=True, exist_ok=True)
    examples = []
    for i, file in enumerate(fa_files, start=1):
        examples.append(f"{example_name}{i}")
    example_number = len(examples)
    ## genome_to_id
    id_to_genome_file = sim_dir_path / "genome_to_id.tsv"
    if not id_to_genome_file.exists():
        with open(id_to_genome_file, "w") as f:
            for i, example in enumerate(examples):
                f.write(f"{example}\t{fa_files[i]}\n")

    ## metadata
    metadata_file = sim_dir_path / "metadata.tsv"
    if not metadata_file.exists():
        with open(metadata_file, "w") as f:
            f.write("genome_ID\tOTU\tNCBI_ID\tnovelty_category\n")
            for i, example in enumerate(examples):
                f.write(f"{example}\t562\t562\tknown_strain\n")   

    ## distribution
    distribution_file = sim_dir_path / "distrubution.txt"
    if not distribution_file.exists():
        with open(distribution_file, "w") as f:
            for i, example in enumerate(examples): 
                f.write(f"{example}\t{distribution[i]}\n")

    ## config_ini
    out_config = sim_dir_path / "config.ini"
    strings = []
    if not out_config.exists():    
        with open(config_ini, "r") as f_in, open(out_config, "w") as f_out:
            for line in f_in:
                if line.strip().startswith("id_to_genome_file"):
                    string = f"id_to_genome_file={id_to_genome_file}"
                    strings.append(string)
                elif line.strip().startswith("distribution_file_paths"):
                    string = f"distribution_file_paths={distribution_file}"
                    strings.append(string)
                elif line.strip().startswith("metadata"):
                    string = f"metadata={metadata_file}"
                    strings.append(string)
                elif line.strip().startswith("size"):
                    string = f"size={size}"
                    strings.append(string)
                elif line.strip().startswith("genomes_total"):
                    string = f"genomes_total={example_number}"
                    strings.append(string)
                elif line.strip().startswith("num_real_genomes"):
                    string = f"num_real_genomes={example_number}"
                    strings.append(string)
                elif line.strip().startswith("max_strains_per_otu"):
                    string = f"max_strains_per_otu={example_number}"
                    strings.append(string)

                elif line.strip().startswith("output_directory"):
                    tokens = line.strip().split("=")
                    sample_id = tokens[1]
                    sim_out_dir = sim_dir_path / tokens[1]
                    if sim_out_dir.exists():
                        raise ValueError(f"The {sim_out_dir} exists, maybe finish simulation.")
                    strings.append(f"output_directory={sim_out_dir}")
                else:
                    strings.append(line.strip())
            f_out.write("\n".join(strings) + "\n") 


rule get_camisim:
    output:
        f"{outdir}/sim_reads/config.ini"
    params:
        distribution = config["distribution"],
        config_ini = config["config_ini"],
        size = config["size"],
    run:
        create_camisim_config(f"{outdir}/sim_reads", fa_files, params.distribution, params.config_ini, params.size)

rule run_camisim:
    input:
        config_ini=rules.get_camisim.output
    params:
        camisim = config["camisim_path"]
    conda:
        config["camisim_env"]
    log:
        camisim_log=f"{outdir}/sim_reads/camisim.log",
        time_log=f"{outdir}/sim_reads/camisim_time.log",
    shell:
        """
        /usr/bin/time -v -o {log.time_log} python {params.camisim} {input.config_ini} > {log.camisim_log} 2>&1
        """ 

rule align:
    input:
        ref_fa=rules.copy_ref.output.ref_fa,
    output: 
        bam = f"{outdir}/align/{example_name}_align.bam",
        bai = f"{outdir}/align/{example_name}_align.bam.bai"
    log:
        "logs/minimap2_align.log"
    params:
        align_threads = config["align_threads"],
        work_threads = config["work_threads"]
    conda:
        vc_env
    shell: 
        """
        mkdir -p {outdir}/align
        read_path=$(ls {outdir}/sim_reads/*/*/reads/*fq.gz)
        minimap2 -ax map-ont {input.ref_fa} $read_path -t {params.align_threads} 2>{log} | \
         samtools view -hb | \
         samtools sort -T {example_name}_align -@ {params.work_threads} > {output.bam} 
        samtools index {output.bam}
        """

rule freebayes_vc:
    input: 
        bam = rules.align.output.bam,
        ref_fa=rules.copy_ref.output.ref_fa,
        fa_fai = f"{outdir}/reference/ref.fa.fai",
    output: 
        vcf = f"{outdir}/tools_result/freebayes/{example_name}.vcf.gz"
    params:
        align_threads = config["align_threads"]
    conda:
        vc_env
    log:
        vc_log="logs/freebayes_vc.log",
        vc_time_log="logs/freebayes_vc_time.log"
    shell:
        """
        mkdir -p {outdir}/tools_result/freebayes
        fasta_generate_regions.py {input.fa_fai} 100000 > {outdir}/tools_result/freebayes/regions.txt
        /usr/bin/time -v -o {log.vc_time_log} freebayes-parallel {outdir}/tools_result/freebayes/regions.txt {params.align_threads} -f {input.ref_fa} {input.bam} --pooled-continuous | bgzip > {output.vcf} 2>{log.vc_log}
        """ 


rule gatk_vc:
    input: 
        bam = rules.align.output.bam,
        ref_fa=rules.copy_ref.output.ref_fa,
    output: 
        vcf = f"{outdir}/tools_result/gatk/{example_name}.vcf.gz"
    params:
        align_threads = config["align_threads"]
    conda:
        vc_env
    log:
        vc_log="logs/gatk_vc.log",
        vc_time_log="logs/gatk_vc_time.log"
    shell:
        """
        mkdir -p {outdir}/tools_result/gatk
        samtools addreplacerg -r "@RG\tID:1\tLB:lib1\tPL:ONT\tPU:unit1\tSM:{example_name}" \
            -o {outdir}/tools_result/gatk/{example_name}_withRG.bam {input.bam}
        samtools index {outdir}/tools_result/gatk/{example_name}_withRG.bam
        /usr/bin/time -v -o {log.vc_time_log} gatk Mutect2 \
            -R {input.ref_fa} \
            -I {outdir}/tools_result/gatk/{example_name}_withRG.bam \
            -O {output.vcf} \
            --native-pair-hmm-threads {params.align_threads} 2>{log.vc_log}       
        """

rule longshot_vc:
    input: 
        bam = rules.align.output.bam,
        ref_fa=rules.copy_ref.output.ref_fa,
    output: 
        vcf = f"{outdir}/tools_result/longshot/{example_name}.vcf.gz"
    conda:
        vc_env
    log:
        vc_log="logs/longshot_vc.log",
        vc_time_log="logs/longshot_vc_time.log"
    shell:
        """
        mkdir -p {outdir}/tools_result/longshot
        /usr/bin/time -v -o {log.vc_time_log} longshot --bam {input.bam} --ref {input.ref_fa} --out {outdir}/tools_result/longshot/{example_name}.vcf --out_bam {outdir}/tools_result/longshot/{example_name}.bam -c 2 -n 2>{log.vc_log}
        bgzip {outdir}/tools_result/longshot/{example_name}.vcf
        """

rule snippy_vc:
    input: 
        bam = rules.align.output.bam,
        ref_fa=rules.copy_ref.output.ref_fa
    output: 
        vcf = f"{outdir}/tools_result/snippy/{example_name}.vcf.gz"
    conda:
        vc_env2
    log:
        vc_log="logs/snippy_vc.log",
        vc_time_log="logs/snippy_vc_time.log"
    shell:
        """
        mkdir -p {outdir}/tools_result/snippy
        samtools addreplacerg -r "@RG\tID:1\tLB:lib1\tPL:ONT\tPU:unit1\tSM:{example_name}" \
            -o {outdir}/tools_result/snippy/{example_name}_withRG.bam {input.bam}
        /usr/bin/time -v -o {log.vc_time_log} snippy --outdir {outdir}/tools_result/snippy/{example_name} --ref {input.ref_fa} --bam {outdir}/tools_result/snippy/{example_name}_withRG.bam 2>{log.vc_log} 
        mv {outdir}/tools_result/snippy/{example_name}/snps.vcf.gz {output.vcf}
        """

rule clair3_vc:
    input: 
        bam = rules.align.output.bam,
        ref_fa=rules.copy_ref.output.ref_fa
    output: 
        vcf = f"{outdir}/tools_result/clair3/{example_name}.vcf.gz"
    params:
        model = config["clair3_ontR10_model"],
        align_threads = config["align_threads"]
    conda:
        vc_env3
    log:
        vc_log="logs/clair3_vc.log",
        vc_time_log="logs/clair3_vc_time.log"
    shell:
        """
        mkdir -p {outdir}/tools_result/clair3
        /usr/bin/time -v -o {log.vc_time_log} run_clair3.sh \
          --bam_fn={input.bam} \
          --ref_fn={input.ref_fa} \
          --threads={params.align_threads} \
          --platform="ont" \
          --model_path={params.model} \
          --output={outdir}/tools_result/clair3/clair_res \
          --include_all_ctgs >{log.vc_log} 2>&1     
        cp -f {outdir}/tools_result/clair3/clair_res/full_alignment.vcf.gz {output.vcf}
        """

rule eval:
    input: 
        vcf = outdir + "/tools_result/{vc_tool}/" + example_name + ".vcf.gz",
        truth = rules.process_vcf.output.final_vcf,
        sdf = f"{outdir}/reference/ref.sdf"
    output: 
        eval_txt = outdir + "/tools_result/{vc_tool}/rtg_eval.txt"
    conda:
        vc_env
    shell:
        """
        ## metagenome should not consider genotypes
        if [ ! -f {input.vcf}.tbi ];then
            tabix -p vcf {input.vcf}
        fi
        if [ -d {outdir}/tools_result/{wildcards.vc_tool}/benchmark ];then
            rm -rf {outdir}/tools_result/{wildcards.vc_tool}/benchmark
        fi
        rtg vcfeval \
            -c {input.vcf} \
            -b {input.truth} \
            -o {outdir}/tools_result/{wildcards.vc_tool}/benchmark \
            -t {input.sdf} \
            --squash-ploidy 1> {output.eval_txt}
        """

rule eval_alt2:
    input: 
        vcf = outdir + "/tools_result/{vc_tool}/" + example_name + ".vcf.gz",
        truth = rules.filter_alt_counts.output.alt2_vcf,
        sdf = f"{outdir}/reference/ref.sdf"
    output: 
        eval_txt = outdir + "/tools_result/{vc_tool}/rtg_eval_alt2.txt"
    conda:
        vc_env
    shell:
        """
        ## metagenome should not consider genotypes
        if [ ! -f {input.vcf}.tbi ];then
            tabix -p vcf {input.vcf}
        fi
        if [ -d {outdir}/tools_result/{wildcards.vc_tool}/benchmark_alt2 ];then
            rm -rf {outdir}/tools_result/{wildcards.vc_tool}/benchmark_alt2
        fi
        rtg vcfeval \
            -c {input.vcf} \
            -b {input.truth} \
            -o {outdir}/tools_result/{wildcards.vc_tool}/benchmark_alt2 \
            -t {input.sdf} \
            --squash-ploidy 1> {output.eval_txt}
        """

rule eval_alt3:
    input: 
        vcf = outdir + "/tools_result/{vc_tool}/" + example_name + ".vcf.gz",
        truth = rules.filter_alt_counts.output.alt3_vcf,
        sdf = f"{outdir}/reference/ref.sdf"
    output: 
        eval_txt = outdir + "/tools_result/{vc_tool}/rtg_eval_alt3.txt"
    conda:
        vc_env
    shell:
        """
        ## metagenome should not consider genotypes
        if [ ! -f {input.vcf}.tbi ];then
            tabix -p vcf {input.vcf}
        fi
        if [ -d {outdir}/tools_result/{wildcards.vc_tool}/benchmark_alt3 ];then
            rm -rf {outdir}/tools_result/{wildcards.vc_tool}/benchmark_alt3
        fi
        rtg vcfeval \
            -c {input.vcf} \
            -b {input.truth} \
            -o {outdir}/tools_result/{wildcards.vc_tool}/benchmark_alt3 \
            -t {input.sdf} \
            --squash-ploidy 1> {output.eval_txt}
        """

rule eval2:
    input: 
        vcf = outdir + "/tools_result/{vc_tool}/" + example_name + ".vcf.gz",
        truth = rules.process_vcf.output.final_vcf,
    output: 
        eval_txt = outdir + "/tools_result/{vc_tool}/eval.txt"
    conda:
        vc_env
    # log:
    #     "logs/eval_{vc_tool}.log"
    
    shell:
        """
        python scripts/eval.py -t {input.truth} -i {input.vcf} > {output.eval_txt}
        python scripts/eval.py -t {input.truth} -i {input.vcf} --po >> {output.eval_txt}
        """

rule hap_py_eval:
    input: 
        vcf = outdir + "/tools_result/{vc_tool}/" + example_name + ".vcf.gz",
        truth = rules.process_vcf.output.final_vcf,      
        ref_fa=rules.copy_ref.output.ref_fa
    output: 
        eval_txt = "{outdir}/tools_result/{vc_tool}/hap_py_eval.txt"
    conda:
        "hap_py"
    shell: 
        """
        som.py \
          {input.truth} \
          {input.vcf} \
          -o {outdir}/tools_result/{wildcards.vc_tool}/benchmark_som_py \
          -r {input.ref_fa} >> {output.eval_txt}       
        """

rule rtg_eval_single:
    input: 
        vcf = outdir + "/tools_result/{vc_tool}/" + example_name + ".vcf.gz",
        truth = outdir + "/simulated/{vcf_name_single}.refseq2simseq.SNP.seq.vcf",
        sdf = f"{outdir}/reference/ref.sdf"
    output: 
        eval_txt = outdir + "/tools_result/{vc_tool}/rtg_eval_{vc_tool}_{vcf_name_single}.txt"
    conda:
        vc_env
    shell:
        """
        if [ ! -f {input.truth}.gz ];then
            bgzip {input.truth} > {input.truth}.gz
        fi

        if [ ! -f {input.truth}.gz.tbi ];then
            tabix -p vcf {input.truth}.gz
        fi

        if [ -d {outdir}/tools_result/{wildcards.vc_tool}/single/benchmark_{wildcards.vcf_name_single} ];then
            rm -rf {outdir}/tools_result/{wildcards.vc_tool}/single/benchmark_{wildcards.vcf_name_single}
        fi
        rtg vcfeval \
            -c {input.vcf} \
            -b {input.truth}.gz \
            -o {outdir}/tools_result/{wildcards.vc_tool}/single/benchmark_{wildcards.vcf_name_single} \
            -t {input.sdf} \
            --squash-ploidy --sample ALT 1> {output.eval_txt}        
        """

rule my_eval_single:
    input: 
        vcf = outdir + "/tools_result/{vc_tool}/" + example_name + ".vcf.gz",
        truth = outdir + "/simulated/{vcf_name_single}.refseq2simseq.SNP.seq.vcf",
        sdf = f"{outdir}/reference/ref.sdf"
    output: 
        eval_txt = outdir + "/tools_result/{vc_tool}/my_eval_{vc_tool}_{vcf_name_single}.txt"
    conda:
        vc_env
    shell:
        """
        python scripts/eval.py -t {input.truth} -i {input.vcf} > {output.eval_txt}
        python scripts/eval.py -t {input.truth} -i {input.vcf} --po >> {output.eval_txt}     
        """

def read_rtg_eval(file):
    precision = recall = 0
    with open(file, "r") as f:
        for i, line in enumerate(f):
            if "None" in line:
                tokens = line.strip().split()
                if len(tokens) < 8: print(file)
                precision = tokens[-3]
                recall = tokens[-2]
    return precision, recall

def read_hap_eval(file):
    hap_precision = hap_recall = 0
    with open(file, "r") as f:
        for line in f:
            if "records" in line:
                tokens = line.strip().split()
                hap_recall = tokens[9]
                hap_precision = tokens[13]
    return hap_precision, hap_recall

def read_my_eval(file):
    metrics = []
    with open(file, "r") as f:
        for line in f:
            if "Precision" in line or "Recall" in line:
                m = line.strip().split(":")[1].strip()
                metrics.append(m)
    return metrics

def get_cov(file): 
    coverage_list = []
    record = False
    with open(file, "r") as f:
        for line in f:
            if "Simulation parameters" in line:
                record = True
            # if line.strip().startswith("Fold Coverage"):
            elif line.strip().startswith("depth") and record:
                coverage = line.split(":")[1].strip().split("X")[0]
                coverage_float = float(coverage)
                coverage_int = int(round(coverage_float))
                coverage_list.append(coverage_int)
                record = False
    coverage_list_sorted = sorted(coverage_list)
    coverage_list_sorted = [str(i) for i in coverage_list_sorted]
    return ",".join(coverage_list_sorted)

def write_eval_file(tool, files, i):
    rtg_precision = rtg_recall = hap_precision = hap_recall = my1_precision = my1_recall = my2_precision = my2_recall = alt2_rtg_recall = alt3_rtg_recall = 0
    single_rtg_evals = []
    single_my_evals1 = []
    single_my_evals2 = []
    for file in files:
        if "rtg_eval.txt" in file:
            rtg_precision, rtg_recall = read_rtg_eval(file)
        elif "hap" in file:
            hap_precision, hap_recall = read_hap_eval(file)
        elif "eval.txt" in file:
            metrics = read_my_eval(file)
            my1_precision, my1_recall, my2_precision, my2_recall = metrics[0], metrics[1], metrics[2], metrics[3]
        elif "alt2" in file:
            _, alt2_rtg_recall = read_rtg_eval(file)
        elif "alt3" in file:
            _, alt3_rtg_recall = read_rtg_eval(file)
        elif "seed" in file and "rtg" in file:
            _, single_recall = read_rtg_eval(file)
            single_rtg_evals.append(single_recall)
        elif "seed" in file and "rtg" not in file:
            metrics = read_my_eval(file)
            single_my_evals1.append(metrics[1])
            single_my_evals2.append(metrics[3])

    num_columns = 14
    row = "|" + "|".join([":----:"] * num_columns) + "|"
    single_rtg_evals_list = ",".join(single_rtg_evals)
    single_my_evals1_list = ",".join(single_my_evals1)
    single_my_evals2_list = ",".join(single_my_evals2)

    cov_s = get_cov(f"{outdir}/sim_reads/camisim.log")
    with open(f"{outdir}/all_eval.md", "a") as f:
        if i == 0:
            f.write(f"tool | rtg_P | rtg_R | hap_P | hap_R | my1_P | my1_R | my2_P | my2_R | alt2_rtg_R | alt3_rtg_R | single_rtg ({cov_s}) | single_my1 ({cov_s}) | single_my2 ({cov_s}) |\n")
            f.write(row + "\n")
        # f.write(f"{rtg_precision}\t{rtg_recall}\t{hap_precision}\t{hap_recall}\t{my1_precision}\t{my1_recall}\t{my2_precision}\t{my2_recall}\t{alt2_rtg_recall}\t{alt3_rtg_recall}\t{single_rtg_evals_list}\t{single_my_evals1_list}\t{single_my_evals2_list}\n")
        f.write(f"|{tool}|{rtg_precision}|{rtg_recall}|{hap_precision}|{hap_recall}|{my1_precision}|{my1_recall}|{my2_precision}|{my2_recall}|{alt2_rtg_recall}|{alt3_rtg_recall}|{single_rtg_evals_list}|{single_my_evals1_list}|{single_my_evals2_list}|\n")
        
            

rule print_res:
    input:
        expand(outdir + "/tools_result/{vc_tool}/eval.txt", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/rtg_eval.txt", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/rtg_eval_alt2.txt", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/rtg_eval_alt3.txt", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/hap_py_eval.txt", vc_tool=vc_tools),
        expand(outdir + "/tools_result/{vc_tool}/rtg_eval_{vc_tool}_{vcf_name_single}.txt", vc_tool=vc_tools, vcf_name_single=snp_vcfs_name),
        expand(outdir + "/tools_result/{vc_tool}/my_eval_{vc_tool}_{vcf_name_single}.txt", vc_tool=vc_tools, vcf_name_single=snp_vcfs_name)
    output: 
        f"{outdir}/all_eval.md"
    run:
        from collections import defaultdict
        res_dict = defaultdict(list)
        for file in input:
            tool_name = file.split("/")[-2]
            res_dict[tool_name].append(file)
        for i, (tool_name, files) in enumerate(res_dict.items()):
            print(f"Processing tool: {tool_name}, index: {i}")
            print(f"files: {files}")
            print(f"Number of files: {len(files)}")
            write_eval_file(tool_name, files, i)


