import argparse
from os.path import join
from .readsim import fqsim
from .lite_utils import utilities

class sim_all(object):

    def __init__(self):
        self.parser = argparse.ArgumentParser(prog='downsample',
                                 description='generating fastq using the downsampled csv',
                                 epilog='_______________________________________________'
                                )
        ### Required Parameters
        self.parser.add_argument('--sample', required=True, type=str, help='sample name (director containing downsampled csv in the path)')
        self.parser.add_argument('--thread', required=True, type=int, help='maximum thread to run')
        self.parser.add_argument('--path', required=True, type=str, help='output path of simulated data')
        self.parser.add_argument(
            '--csv',
            required=False,
            type=str,
            help='path to a molecule pool CSV/TSV (default: {path}/{sample}/{sample}.lib.csv).',
        )
        self.parser.add_argument(
            '--seed',
            required=False,
            type=int,
            default=None,
            help='set a random seed to repeat result.(default=None)',
        )

        ### Fastq parameters
        self.parser.add_argument(
            '--skip-hifi',
            required=False,
            action='store_true',
            help='skip PacBio HiFi long-read simulation. (optional)',
        )
        self.parser.add_argument(
            '--skip-ont',
            required=False,
            action='store_true',
            help='skip Oxford Nanopore (ONT) long-read simulation. (optional)',
        )
        self.parser.add_argument('--ont-model', required=False, type=str, default='R94', help='paras for PBSIM2 option(R94, R95, R103, P4C2, P5C3, P6C4)')
        self.parser.add_argument('--ont-mean', required=False, type=float, default=3000, help='mean of reads length of simulated ONT fastq (default=3000)')
        self.parser.add_argument('--ont-std', required=False, type=float, default=2500, help='std  of reads length of simulated ONT fastq (default=2500)')
        self.parser.add_argument(
            '--hifi-sample-fastq',
            required=False,
            type=str,
            help='HiFi sample fastq for PBSIM2 sampling-based simulation. (used as --sample-fastq)',
        )
        self.parser.add_argument(
            '--hifi-mode',
            required=False,
            type=str,
            default='auto',
            choices=['auto', 'sampling', 'simple'],
            help='HiFi simulation mode: auto (sampling if profile/sample is provided, else simple), sampling (PBSIM2 sampling-based), simple (built-in synthetic model). (default=auto)',
        )
        self.parser.add_argument(
            '--hifi-profile-id',
            required=False,
            type=str,
            help='HiFi profile id for PBSIM2 sampling-based simulation. (optional)',
        )
        self.parser.add_argument(
            '--hifi-profile-root',
            required=False,
            type=str,
            help='directory used as HOME for PBSIM2 HiFi sampling profiles (portable cache). (optional)',
        )
        self.parser.add_argument(
            '--hifi-len-min',
            required=False,
            type=int,
            default=5000,
            help='HiFi simple-mode minimum read length. (default=5000)',
        )
        self.parser.add_argument(
            '--hifi-len-peak-min',
            required=False,
            type=int,
            default=10000,
            help='HiFi simple-mode peak-range minimum read length. (default=10000)',
        )
        self.parser.add_argument(
            '--hifi-len-peak-max',
            required=False,
            type=int,
            default=25000,
            help='HiFi simple-mode peak-range maximum read length. (default=25000)',
        )
        self.parser.add_argument(
            '--hifi-len-max',
            required=False,
            type=int,
            default=60000,
            help='HiFi simple-mode maximum read length. (default=60000)',
        )
        self.parser.add_argument(
            '--hifi-qmin',
            required=False,
            type=int,
            default=20,
            help='HiFi simple-mode minimum Phred quality. (default=20)',
        )
        self.parser.add_argument(
            '--hifi-qmean',
            required=False,
            type=int,
            default=30,
            help='HiFi simple-mode mean Phred quality (also used as constant quality when qsd=0). (default=30)',
        )
        self.parser.add_argument(
            '--hifi-qsd',
            required=False,
            type=float,
            default=0.0,
            help='HiFi simple-mode Phred quality standard deviation. (default=0.0)',
        )
        self.parser.add_argument(
            '--hifi-total-reads',
            required=False,
            type=int,
            help='total number of HiFi reads to generate (simple-mode only). (optional)',
        )
        self.parser.add_argument(
            '--skip-sr',
            required=False,
            action='store_true',
            help='skip short-read (Illumina) simulation. (optional)',
        )
        self.parser.add_argument('--sr-platform', required=False, type=str, default='HS25', help='paras for art option(HS10, HS20,HS25, HSXn, HSXt, MinS, MSv1, MSv3, NS50)')
        self.parser.add_argument('--sr-mean', required=False, type=float, default=400, help='mean of insert length of simulated NGS fastq (default=400)')
        self.parser.add_argument('--sr-std', required=False, type=float, default=125, help='std  of insert length of simulated NGS fastq (default=125)')
        self.parser.add_argument('--sr-readlen', required=False, type=float, default=150, help='reads length of simulated NGS fastq (default=125)')
        self.args = self.parser.parse_args()
        if (
            (not self.args.skip_hifi)
            and (self.args.hifi_mode == 'sampling')
            and (not self.args.hifi_sample_fastq)
            and (not self.args.hifi_profile_id)
        ):
            self.parser.error('--hifi-sample-fastq or --hifi-profile-id is required when HiFi is enabled with sampling mode.')

        self.simulate()

    def simulate(self):
        csv_path = self.args.csv if self.args.csv else join(self.args.path, self.args.sample, self.args.sample + '.lib.csv')
        temp = fqsim(sample = self.args.sample,
                     csv = csv_path,
                     path = self.args.path,
                     seed = self.args.seed,
                     thread = self.args.thread,
                     skip_sr = self.args.skip_sr,
                     skip_hifi = self.args.skip_hifi,
                     skip_ont = self.args.skip_ont,
                     ont_model = self.args.ont_model,
                     ont_mean = self.args.ont_mean,
                     ont_std = self.args.ont_std,
                     hifi_sample_fastq = self.args.hifi_sample_fastq,
                     hifi_mode = self.args.hifi_mode,
                     hifi_profile_id = self.args.hifi_profile_id,
                     hifi_profile_root = self.args.hifi_profile_root,
                     hifi_len_min = self.args.hifi_len_min,
                     hifi_len_peak_min = self.args.hifi_len_peak_min,
                     hifi_len_peak_max = self.args.hifi_len_peak_max,
                     hifi_len_max = self.args.hifi_len_max,
                     hifi_qmin = self.args.hifi_qmin,
                     hifi_qmean = self.args.hifi_qmean,
                     hifi_qsd = self.args.hifi_qsd,
                     hifi_total_reads = self.args.hifi_total_reads,
                     sr_platform = self.args.sr_platform,
                     sr_mean = self.args.sr_mean,
                     sr_std = self.args.sr_std,
                     sr_readlen = self.args.sr_readlen,
                     )

def main():
    run = sim_all()

if __name__ == '__main__':
    main()
