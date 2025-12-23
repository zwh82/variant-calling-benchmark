
import sys, gzip, argparse

usage = "variant calling benchmark"

def main():
    parser = argparse.ArgumentParser(prog="python genomes_process.py", description=usage)
    parser.add_argument("-i", "--pred", dest="pred", type=str, help="predicted vcf file.")
    parser.add_argument("-t", "--truth", dest="truth", type=str, help="truth vcf file.")
    parser.add_argument("-b", "--base", dest="baseline", type=str, help="baseline vcf file from rtg.")
    parser.add_argument("--sep", dest="separate", action="store_true", help="separate alts into single and then eval.")
    parser.add_argument("--sort", dest="sort", action="store_true", help="sort alts.")
    parser.add_argument("--po", dest="pos_only", action="store_true", help="eval only considers pos.")
    parser.add_argument("--tool", dest="print_tool", action="store_true", help="print tool.")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    eval_main(args)

def eval_main(args):
    if args.print_tool:
        try:
            tool = pred_vcf.split("/")[2]
            print(tool)
        except:
            pass
        
    if args.separate:
        eval_single(args.truth, args.pred, args.pos_only)
    else:
        eval(args.truth, args.pred, args.baseline, args.pos_only, args.sort)

def eval(truth_vcf, pred_vcf, baseline_vcf, pos_only, is_sort):
    truth = load_variants(truth_vcf, pos_only, is_sort)
    pred = load_variants(pred_vcf, pos_only, is_sort)
    tp, fp, fn, precision, recall = compute_metrics(truth, pred)
    if baseline_vcf:
        baseline, variants_dict = load_variants2(baseline_vcf)
        diff_baseline_tp(tp, baseline, variants_dict)

    if pos_only:
        print("eval only considers pos")
    else:
        print("eval considers pos, ref, alt")

    print(f"Number of truth variants: {len(truth)}")
    print(f"Number of predicted variants: {len(pred)}")
    print(f"True positives (TP): {len(tp)}")
    print(f"False positives (FP): {len(fp)}")
    print(f"False negatives (FN): {len(fn)}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print("-"*50)

def load_variants(vcf_file, pos_only, is_sort):
    variants = set()
    with open_file(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            # alt = parts[4] 
            if is_sort:
                alt_list = parts[4].split(",")
                sorted_alt_list = sorted(alt_list) 
                alt = ",".join(sorted_alt_list)
            else:
                alt = parts[4] 
            if pos_only:
                variants.add((chrom, pos))
            else:
                variants.add((chrom, pos, ref, alt))

    return variants

def load_variants2(vcf_file):
    variants = set()
    variants_dict = dict()
    with open_file(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt = parts[4] 
            # variants.add((chrom, pos, ref, alt))
            variants.add((chrom, pos))
            variants_dict[(chrom, pos)] = (ref, alt)
    return variants, variants_dict

def load_variants3(vcf_file, pos_only):
    """
    all alts is split into a single alt, and then append to variant set
    """
    variants = set()
    with open_file(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alts = parts[4].split(",") 
            for alt in alts:
                if alt == "": continue
                if pos_only:
                    variants.add((chrom, pos))
                else:
                    variants.add((chrom, pos, ref, alt))
    return variants

def eval_single(truth_vcf, pred_vcf, pos_only):
    truth = load_variants3(truth_vcf, pos_only)
    pred = load_variants3(pred_vcf, pos_only)
    recall = compute_metrics3(truth, pred)
    if pos_only:
        print("eval only considers pos")
    else:
        print("eval considers pos, ref, alt") 
    print(f"recall: {recall}\n")

def compute_metrics(truth_set, pred_set):
    tp = truth_set & pred_set
    fp = pred_set - truth_set
    fn = truth_set - pred_set

    precision = len(tp) / (len(tp) + len(fp)) if (len(tp) + len(fp)) > 0 else 0.0
    recall = len(tp) / (len(tp) + len(fn)) if (len(tp) + len(fn)) > 0 else 0.0

    return tp, fp, fn, precision, recall

def compute_metrics3(truth_set, pred_set):
    tp = truth_set & pred_set
    recall = len(tp) / len(truth_set)
    return recall

def diff_baseline_tp(tp, baseline, variants_dict):
    diff1 = baseline - tp
    with open("diff_baseline_tp.txt", "w") as f:
        for _d in diff1:
            ref_alt = variants_dict[_d]
            f.write(f"{_d[0]}\t{_d[1]}\t{ref_alt[0]}\t{ref_alt[1]}\n")   
 

def open_file(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")

if __name__ == "__main__":
    sys.exit(main())
