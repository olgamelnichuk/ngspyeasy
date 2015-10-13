#!/usr/bin/env python
import os


def select(gb, projects_home):
    if gb == "b37":
        return GenomeBuild(
            projects_home.resource_path("reference_genomes_b37"),
            GENOMEINDEX="human_g1k_v37",
            REFFASTA="human_g1k_v37.fasta",
            DBSNP="dbsnp_138.b37.recab",
            KNOWN_INDELS="Mills_and_1000G_gold_standard.indels.b37.vcf",
            KNOWN_SNPS_1000G="1000G_phase1.snps.high_confidence.b37.vcf",
            KNOWN_SNPS_OMNI="1000G_omni2.5.b37.vcf",
            KNOWN_SNPS_b138="dbsnp_138.b37.vcf"
        )
    elif gb == "hg19":
        return GenomeBuild(
            projects_home.resource_path("reference_genomes_hg19"),
            GENOMEINDEX="ucsc.hg19",
            REFFASTA="ucsc.hg19.fasta",
            DBSNP="dbsnp_138.hg19.recab",
            KNOWN_INDELS="Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
            KNOWN_SNPS_1000G="1000G_phase1.snps.high_confidence.hg19.sites.vcf",
            KNOWN_SNPS_OMNI="1000G_omni2.5.hg19.sites.vcf",
            KNOWN_SNPS_b138="bsnp_138.hg19.vcf"
        )

    elif gb == "hs37d5":
        return GenomeBuild(
            projects_home.resource_path("reference_genomes_hs37d5"),
            GENOMEINDEX="hs37d5"
        )
    elif gb == "hs38DH":
        return GenomeBuild(
            projects_home.resource_path("reference_genomes_hs38DH"),
            GENOMEINDEX="hs38DH"
        )


class GenomeBuild(object):
    def __init__(self, ref_dir, **kwargs):
        self.ref_dir = ref_dir
        self.args = dict().update(kwargs)

    def refdir(self):
        return self.ref_dir

    def path_to(self, value):
        return None if value is None else os.path.join(self.ref_dir, value)

    def genome_index(self):
        return self.path_to(self.args.get("GENOMEINDEX", None))

    def ref_fasta(self):
        return self.path_to(self.args.get("REFFASTA", None))

    def dbsnp(self):
        return self.path_to(self.args.get("DBSNP", None))

    def known_indels(self):
        return self.path_to(self.args.get("KNOWN_INDELS", None))

    def known_snps_100g(self):
        return self.path_to(self.args.get("KNOWN_SNPS_1000G", None))

    def known_snps_omni(self):
        return self.path_to(self.args.get("KNOWN_SNPS_OMNI", None))

    def known_snps_b138(self):
        return self.path_to(self.args.get("KNOWN_SNPS_b138", None))

    def adapter_fa(self):
        return self.path_to("contaminant_list.fa")
