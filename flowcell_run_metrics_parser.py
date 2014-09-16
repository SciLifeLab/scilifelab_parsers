import os
import re
import yaml
import glob
import xml.parsers.expat
#from uuid import uuid4
import json
import numpy as np
import csv
import collections
#import xml.etree.cElementTree as ET
from bs4 import BeautifulSoup
import datetime

#from scilifelab.log import minimal_logger
#LOG = minimal_logger("bcbio")

#from bcbio.broad.metrics import PicardMetricsParser
#from bcbio.pipeline.qcsummary import FastQCParser

class FlowcellRunMetricsParser(RunMetricsParser):
    """Flowcell level class for parsing flowcell run metrics data."""
    _lanes = range(1,9)
    def __init__(self, path):
        RunMetricsParser.__init__(self)
        self.path = path
        self._collect_files()

    def parseRunInfo(self, fn="RunInfo.xml", **kw):
        infile = os.path.join(os.path.abspath(self.path), fn)
        self.log.debug("parseRunInfo: going to read {}".format(infile))
        if not os.path.exists(infile):
            self.log.warn("No such file {}".format(infile))
            return {}
        try:
            fp = open(infile)
            parser = RunInfoParser()
            data = parser.parse(fp)
            fp.close()
            return data
        except:
            self.log.warn("Reading file {} failed".format(os.path.join(os.path.abspath(self.path), fn)))
            return {}

    def parseRunParameters(self, fn="runParameters.xml", **kw):
        """Parse runParameters.xml from an Illumina run.

        :param fn: filename
        :param **kw: keyword argument

        :returns: parsed data structure
        """
        infile = os.path.join(os.path.abspath(self.path), fn)
        self.log.debug("parseRunParameters: going to read {}".format(infile))
        if not os.path.exists(infile):
            self.log.warn("No such files {}".format(infile))
            return {}
        try:
            with open(infile) as fh:
                parser = RunParametersParser()
                data = parser.parse(fh)
            return data
        except:
            self.log.warn("Reading file {} failed".format(os.path.join(os.path.abspath(self.path), fn)))
            return {}

    def parseDemultiplexConfig(self, fn="DemultiplexConfig.xml", **kw):
        """Parse the DemultiplexConfig.xml configuration files"""
        pattern = os.path.join(os.path.abspath(self.path), "Unaligned*", fn)
        cfg = {}
        for cfgfile in glob.glob(pattern):
            parser = DemultiplexConfigParser(cfgfile)
            data = parser.parse()
            if len(data) > 0:
                cfg[os.path.basename(os.path.dirname(cfgfile))] = data
        return cfg

    def parse_samplesheet_csv(self, runinfo_csv="SampleSheet.csv", **kw):
        infile = os.path.join(os.path.abspath(self.path), runinfo_csv)
        self.log.debug("parse_samplesheet_csv: going to read {}".format(infile))
        if not os.path.exists(infile):
            self.log.warn("No such file {}".format(infile))
            return {}
        try:
            fp = open(infile)
            runinfo = [x for x in csv.DictReader(fp)]
            fp.close()
            return runinfo
        except:
            self.log.warn("Reading file {} failed".format(infile))
            return {}

    def parse_run_info_yaml(self, run_info_yaml="run_info.yaml", **kw):
        infile = os.path.join(os.path.abspath(self.path), run_info_yaml)
        self.log.debug("parse_run_info_yaml: going to read {}".format(infile))
        if not os.path.exists(infile):
            self.log.warn("No such file {}".format(infile))
            return {}
        try:
            fp = open(infile)
            runinfo = yaml.load(fp)
            fp.close()
            return runinfo
            return True
        except:
            self.log.warn("No such file {}".format(infile))
            return False

    def parse_illumina_metrics(self, fullRTA=False, **kw):
        self.log.debug("parse_illumina_metrics")
        fn = []
        for root, dirs, files in os.walk(os.path.abspath(self.path)):
            for f in files:
                if f.endswith(".xml"):
                    fn.append(os.path.join(root, f))
        self.log.debug("Found {} RTA files {}...".format(len(fn), ",".join(fn[0:10])))
        parser = IlluminaXMLParser()
        metrics = parser.parse(fn, fullRTA)
        def filter_function(f):
            return f is not None and f == "run_summary.json"
        try:
            metrics.update(self.parse_json_files(filter_fn=filter_function).pop(0))
        except IndexError:
            pass
        return metrics

    def parse_filter_metrics(self, fc_name, **kw):
        """pre-CASAVA: Parse filter metrics at flowcell level"""
        self.log.debug("parse_filter_metrics for flowcell {}".format(fc_name))
        lanes = {str(k):{} for k in self._lanes}
        for lane in self._lanes:
            pattern = "{}_[0-9]+_[0-9A-Za-z]+(_nophix)?.filter_metrics".format(lane)
            lanes[str(lane)]["filter_metrics"] = {"reads":None, "reads_aligned":None, "reads_fail_align":None}
            files = self.filter_files(pattern)
            self.log.debug("filter metrics files {}".format(",".join(files)))
            try:
                fp = open(files[0])
                parser = MetricsParser()
                data = parser.parse_filter_metrics(fp)
                fp.close()
                lanes[str(lane)]["filter_metrics"] = data
            except:
                self.log.warn("No filter nophix metrics for lane {}".format(lane))
        return lanes

    def parse_bc_metrics(self, fc_name, **kw):
        """Parse bc metrics at sample level"""
        self.log.debug("parse_bc_metrics for flowcell {}".format(fc_name))
        lanes = {str(k):{} for k in self._lanes}
        for lane in self._lanes:
            pattern = "{}_[0-9]+_[0-9A-Za-z]+(_nophix)?[\._]bc[\._]metrics".format(lane)
            lanes[str(lane)]["bc_metrics"] = {}
            files = self.filter_files(pattern)
            self.log.debug("bc metrics files {}".format(",".join(files)))
            try:
                parser = MetricsParser()
                fp = open(files[0])
                data = parser.parse_bc_metrics(fp)
                fp.close()
                lanes[str(lane)]["bc_metrics"] = data
            except:
                self.log.warn("No bc_metrics info for lane {}".format(lane))
        return lanes

    def parse_undemultiplexed_barcode_metrics(self, fc_name, **kw):
        """Parse the undetermined indices top barcodes materics
        """

        # Use a glob to allow for multiple fastq folders
        metrics_file_pattern = os.path.join(self.path, "Unaligned*", 
                                    "Basecall_Stats_*{}".format(fc_name[1:]), 
                                    "Undemultiplexed_stats.metrics")
        metrics = {'undemultiplexed_barcodes': []}
        for metrics_file in glob.glob(metrics_file_pattern):
            self.log.debug("parsing {}".format(metrics_file))
            if not os.path.exists(metrics_file):
                self.log.warn("No such file {}".format(metrics_file))
                continue

            with open(metrics_file) as fh:
                parser = MetricsParser()
                in_handle = csv.DictReader(fh, dialect=csv.excel_tab)
                data = parser.parse_undemultiplexed_barcode_metrics(in_handle)
                for lane, items in data.items():
                    for item in items:
                        item['lane'] = lane
                        metrics['undemultiplexed_barcodes'].append(item)

        # Define a function for sorting values according to lane and yield
        def by_lane_yield(data):
            return '{}-{}'.format(data.get('lane',''),data.get('count','').zfill(10))

        # Remove duplicate entries resulting from multiple stats files
        for metric in ['undemultiplexed_barcodes']:
            dedupped = {}
            for row in metrics[metric]:
                key = "\t".join(row.values())
                if key not in dedupped:
                    dedupped[key] = row
                else:
                    self.log.warn("Duplicates of Undemultiplexed barcode entries"
                                " discarded: {}".format(key[0:min(35,len(key))]))

            # Reformat the structure of the data to fit the downstream processing
            lanes = {}
            for row in sorted(dedupped.values(), key=by_lane_yield, reverse=True):
                lane = row['lane']
                if lane not in lanes:
                    lanes[lane] = {metric: {k:[] for k in row.keys()}}
                for k in row.keys():
                    lanes[lane][metric][k].append(row[k])

        return lanes

    def parse_demultiplex_stats_htm(self, fc_name, **kw):
        """Parse the Unaligned*/Basecall_Stats_*/Demultiplex_Stats.htm file
        generated from CASAVA demultiplexing and returns barcode metrics.
        """
        metrics = {"Barcode_lane_statistics": [],
                   "Sample_information": []}
        # Use a glob to allow for multiple fastq directories
        htm_file_pattern = os.path.join(self.path, "Unaligned*", 
                                "Basecall_Stats_*{}".format(fc_name[1:]), 
                                "Demultiplex_Stats.htm")
        for htm_file in glob.glob(htm_file_pattern):
            self.log.debug("parsing {}".format(htm_file))
            if not os.path.exists(htm_file):
                self.log.warn("No such file {}".format(htm_file))
                continue
            with open(htm_file) as fh:
                htm_doc = fh.read()
            soup = BeautifulSoup(htm_doc)
            ##
            ## Find headers
            allrows = soup.findAll("tr")
            column_gen=(row.findAll("th") for row in allrows)
            parse_row = lambda row: row
            headers = [h for h in map(parse_row, column_gen) if h]
            bc_header = [str(x.string) for x in headers[0]]
            smp_header = [str(x.string) for x in headers[1]]
            ## 'Known' headers from a Demultiplex_Stats.htm document
            bc_header_known = ['Lane', 'Sample ID', 'Sample Ref', 'Index', 
                'Description', 'Control', 'Project', 'Yield (Mbases)', '% PF', 
                '# Reads', '% of raw clusters per lane', 
                '% Perfect Index Reads', '% One Mismatch Reads (Index)', 
                '% of >= Q30 Bases (PF)', 'Mean Quality Score (PF)']
            smp_header_known = ['None', 'Recipe', 'Operator', 'Directory']
            if not bc_header == bc_header_known:
                self.log.warn("Barcode lane statistics header information has"
                        " changed. New format?\nOld format: {}\nSaw: {}".format(
                        ",".join((["'{}'".format(x) for x in bc_header_known])),
                        ",".join(["'{}'".format(x) for x in bc_header])))
            if not smp_header == smp_header_known:
                self.log.warn("Sample header information has changed. New "
                        "format?\nOld format: {}\nSaw: {}".format(",".join((
                        ["'{}'".format(x) for x in smp_header_known])), 
                        ",".join(["'{}'".format(x) for x in smp_header])))
            ## Fix first header name in smp_header since htm document is 
            ## mal-formatted: <th>Sample<p></p>ID</th>
            smp_header[0] = "Sample ID"

            ## Parse Barcode lane statistics
            soup = BeautifulSoup(htm_doc)
            table = soup.findAll("table")[1]
            rows = table.findAll("tr")
            column_gen = (row.findAll("td") for row in rows)
            for i in range(0, len(bc_header)) if row}:
                parse_row = lambda row: {bc_header[i]:str(row[i].string) 
            metrics["Barcode_lane_statistics"].extend(map(parse_row, column_gen))

            ## Parse Sample information
            soup = BeautifulSoup(htm_doc)
            table = soup.findAll("table")[3]
            rows = table.findAll("tr")
            column_gen = (row.findAll("td") for row in rows)
            for i in range(0, len(smp_header)) if row}:
                parse_row = lambda row: {smp_header[i]:str(row[i].string)
            metrics["Sample_information"].extend(map(parse_row, column_gen))

        # Define a function for sorting the values
        def by_lane_sample(data):
            return "{}-{}-{}".format(data.get('Lane',''), 
                                data.get('Sample ID',''),data.get('Index',''))

        # Post-process the metrics data to eliminate duplicates resulting from 
        # multiple stats files
        for metric in ['Barcode_lane_statistics', 'Sample_information']:
            dedupped = {}
            for row in metrics[metric]:
                key = "\t".join(row.values())
                if key not in dedupped:
                    dedupped[key] = row
                else:
                    self.log.debug("Duplicates of Demultiplex Stats entries "
                                "discarded: {}".format(key[0:min(35,len(key))]))
            metrics[metric] = sorted(dedupped.values(), key=by_lane_sample)

        ## Set data
        return metrics
