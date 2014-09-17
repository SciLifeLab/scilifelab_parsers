"""bcbio qc module. Parsers for collecting qc metrics."""
import os
import re
import yaml
import glob
import xml.parsers.expat
from uuid import uuid4
import json
import numpy as np
import csv
import collections
import xml.etree.cElementTree as ET
from bs4 import BeautifulSoup
import datetime
import logging as LOG


class MetricsParser():
    """Basic class for parsing metrics"""
    def __init__(self, log=None):
        self.log = LOG
        if log:
            self.log = log

    def parse_bc_metrics(self, in_handle):
        data = {}
        while 1:
            line = in_handle.readline()
            if not line:
                break
            vals = line.rstrip("\t\n\r").split("\t")
            data[vals[0]] = int(vals[1])
        return data

    def parse_filter_metrics(self, in_handle):
        data = {}
        data["reads"] = int(in_handle.readline().rstrip("\n").split(" ")[-1])
        data["reads_aligned"] = int(in_handle.readline().split(" ")[-2])
        data["reads_fail_align"] = int(in_handle.readline().split(" ")[-2])
        return data

    def parse_fastq_screen_metrics(self, in_handle):
        in_handle.readline()
        data = {}
        while 1:
            line = in_handle.readline()
            if not line:
                break
            vals = line.rstrip("\t\n").split("\t")
            data[vals[0]] = {}
            data[vals[0]]["Unmapped"] = float(vals[1])
            data[vals[0]]["Mapped_One_Library"] = float(vals[2])
            data[vals[0]]["Mapped_Multiple_Libraries"] = float(vals[3])
        return data

    def parse_undemultiplexed_barcode_metrics(self, in_handle):

        data = collections.defaultdict(list)
        for line in in_handle:
            data[line['lane']].append(dict([(c,[line[c],''][line[c] is None]) for c in in_handle.fieldnames if c != 'lane']))
        return data

    def parse_bcbb_checkpoints(self, in_handle):

        TIMEFORMAT = "%Y-%m-%dT%H:%M:%S.%f"
        timestamp = []
        for line in in_handle:
            try:
                ts = "{0}Z".format(datetime.datetime.strptime(line.strip(), TIMEFORMAT).isoformat())
                timestamp.append(ts)
            except ValueError:
                pass

        return timestamp

    def parse_software_versions(self, in_handle):
        sver = {}
        for line in in_handle:
            try:
                s = line.split()
                if len(s) == 2:
                    sver[s[0]] = s[1]
            except:
                pass
        return sver


class FlowcellRunMetricsParser():
    """Flowcell level class for parsing flowcell run metrics data."""
    def __init__(self):
        self.log = LOG

    def parse_undemultiplexed_barcode_metrics(self, metrics_file, **kw):
        """Parse the undetermined indices top barcodes materics
        """
        metric = 'undemultiplexed_barcodes'
        metrics = {metric: []}
        with open(metrics_file) as fh:
            parser = MetricsParser()
            in_handle = csv.DictReader(fh, dialect=csv.excel_tab)
            data = parser.parse_undemultiplexed_barcode_metrics(in_handle)
            for lane, items in data.items():
                for item in items:
                    item['lane'] = lane
                    metrics[metric].append(item)

        # Define a function for sorting values according to lane and yield
        def by_lane_yield(data):
            return '{0}-{1}'.format(data.get('lane',''),data.get('count','').zfill(10))

        # Remove duplicate entries resulting from multiple stats files
        dedupped = {}
        for row in metrics[metric]:
            key = "\t".join(row.values())
            if key not in dedupped:
                dedupped[key] = row
            else:
                self.log.warn("Duplicates of Undemultiplexed barcode entries discarded: {0}".format(key[0:min(35,len(key))]))

        # Reformat the structure of the data to fit the downstream processing
        lanes = {}
        for row in sorted(dedupped.values(), key=by_lane_yield, reverse=True):
            lane = row['lane']
            if lane not in lanes:
                lanes[lane] = {metric: dict([(k,[]) for k in row.keys()])}
            for k in row.keys():
                lanes[lane][metric][k].append(row[k])

        return lanes

    def parse_demultiplex_stats_htm(self, htm_file, **kw):
        """Parse the Demultiplex_Stats.htm file
        generated from CASAVA demultiplexing and returns barcode metrics.
        """
        metrics = {"Barcode_lane_statistics": [], "Sample_information": []}
        self.log.debug("parsing {0}".format(htm_file))
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
                'Description', 'Control', 'Project', 'Yield (Mbases)', 
                '% PF', '# Reads', '% of raw clusters per lane', 
                '% Perfect Index Reads', '% One Mismatch Reads (Index)', 
                '% of >= Q30 Bases (PF)', 'Mean Quality Score (PF)']
        smp_header_known = ['None', 'Recipe', 'Operator', 'Directory']
        if not bc_header == bc_header_known:
            self.log.warn("Barcode lane statistics header information has"
                    " changed. New format?\nOld format: {0}\nSaw: {1}".format(
                    ",".join((["'{0}'".format(x) for x in bc_header_known])),
                    ",".join(["'{0}'".format(x) for x in bc_header])))
        if not smp_header == smp_header_known:
            self.log.warn("Sample header information has changed. New "
                    "format?\nOld format: {0}\nSaw: {1}".format(
                    ",".join((["'{0}'".format(x) for x in smp_header_known])),
                    ",".join(["'{0}'".format(x) for x in smp_header])))
        ## Fix first header name in smp_header since htm document is mal-formatted: <th>Sample<p></p>ID</th>
        smp_header[0] = "Sample ID"

        ## Parse Barcode lane statistics
        soup = BeautifulSoup(htm_doc)
        table = soup.findAll("table")[1]
        rows = table.findAll("tr")
        column_gen = (row.findAll("td") for row in rows)
        parse_row = lambda row: dict([(bc_header[i],str(row[i].string)) for i in range(0, len(bc_header)) if row])
        metrics["Barcode_lane_statistics"].extend(map(parse_row, column_gen))

        ## Parse Sample information
        soup = BeautifulSoup(htm_doc)
        table = soup.findAll("table")[3]
        rows = table.findAll("tr")
        column_gen = (row.findAll("td") for row in rows)
        parse_row = lambda row: dict([(smp_header[i],str(row[i].string)) for i in range(0, len(smp_header)) if row])
        metrics["Sample_information"].extend(map(parse_row, column_gen))

        # Define a function for sorting the values
        def by_lane_sample(data):
            return "{0}-{1}-{2}".format(data.get('Lane',''),data.get('Sample ID',''),data.get('Index',''))

        # Post-process the metrics data to eliminate duplicates resulting from multiple stats files
        for metric in ['Barcode_lane_statistics', 'Sample_information']:
            dedupped = {}
            for row in metrics[metric]:
                key = "\t".join(row.values())
                if key not in dedupped:
                    dedupped[key] = row
                else:
                    self.log.debug("Duplicates of Demultiplex Stats entries discarded: {0}".format(key[0:min(35,len(key))]))
            metrics[metric] = sorted(dedupped.values(), key=by_lane_sample)

        ## Set data
        return metrics
