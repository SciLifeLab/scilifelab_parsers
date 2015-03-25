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
        self.flowcell = None

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

    def parse_demultiplex_stats_htm(self, html_file):
        metrics = self._html_tables_to_lists_of_tuples(html_file)    
        new_metrics = {}
        for header, table in metrics.items():
            print header
            new_metrics[header] = []
            for sample in table:
                new_metrics[header].append(dict(sample))
        print new_metrics
        return new_metrics

    def _html_tables_to_lists_of_tuples(self, html):
        html = open(html,'r')
        bs = BeautifulSoup(html)
        metrics = {}
        keys = []
        header = None
        for table in bs.findAll('table'):
            header = table.findPrevious('h2')
            if header:
                header = header.text.replace(' ','_')
            if header not in metrics.keys():
                metrics[header] = []
            if not metrics[header]: 
                for row in table.findChildren('tr'):
                    if row.findChildren('th'):
                        keys = map(lambda v: v.text, row.findChildren('th')) 
                    values = map(lambda v: v.text, row.findChildren('td'))
                    if len(values)==len(keys):
                        metrics[header].append(zip(keys, values))
        return metrics


#html= open('laneBarcode.html','r')
#l=html_tables_to_dicts(html)

#html= open('Demultiplex_Stats.htm','r')
#D=html_tables_to_dicts(html)
