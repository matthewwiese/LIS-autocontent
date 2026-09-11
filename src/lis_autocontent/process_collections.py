"""Populate collections and resources for JBrowse2 and BLAST from a remote Datastore."""

#!/usr/bin/env python3

import os
import sys
import json
import pathlib
import subprocess
from html.parser import HTMLParser
import requests
import yaml

from .datastore_files import NODE_README_FIELDS


class ProcessCollections:
    """Parses Collections from the datastore_url provided. Default: https://data.legumeinfo.org"""

    def __init__(
        self,
        logger=None,
        datastore_url="https://data.legumeinfo.org",
        jbrowse_url="",
        out_dir="./autocontent",
    ):
        self.logger = logger
        if self.logger:
            self.logger.info("logger initialized")
        else:  # logger object required
            print("logger required to initialize process_collections")
            sys.exit(1)
        self.from_github = None  # read from github
        self.collections = []  # stores all collections from self.parse_attributes
        self.datastore_url = datastore_url  # URL to search for collections
        self.jbrowse_url = jbrowse_url  # URL to append jbrowse2 sessions
        self.out_dir = out_dir  # output directory for objects.  This is set by the runtimes if provided
        self.files = (
            {}
        )  # stores all files by collection type. This is used to populate output after scanning
        self.file_objects = []  # a list of all file objects to write for DSCensor nodes
        self.collection_types = [  # collection types currently recorded from datastore_url
            "genomes",  # fasta
            "annotations",  # gff3
            "diversity",  # vcf
            "expression",  # bed, wig, bw
            "genetic",  # bed
            "markers",  # bed
            "synteny",  # paf, this may be depricated now for genome_alignments
            "genome_alignments",  # paf
        ]  # types to search the datastore_url for
        #        self.relationships = {'genomes': {'annotations': ...}, 'annotations': {}}  # establish related objects once this is relevant
        self.current_taxon = {}
        self.species_descriptions = (
            []
        )  # list of all species descriptions to be written to species collections
        self.infraspecies_resources = {}  # used to track all "strains"
        self.species_collections_handle = (
            None  # yaml file to write for species collections
        )
        self.genus_resources_handle = None  # yaml file to write for genus resources
        self.species_resources_handle = None  # yaml file to write for species resources

    def get_remote(self, url):
        """Uses requests.get to grab remote URL returns response.text otherwise returns False"""
        logger = self.logger
        response = requests.get(url, timeout=5)  # get remote object
        if response.status_code == 200:  # SUCCESS
            return response.text
        logger.debug(f"GET failed with status {response.status_code} for: {url}")
        return False

    def head_remote(self, url):
        """Uses requests.head to grab remote URL returns response.text otherwise returns False"""
        logger = self.logger
        response = requests.head(url, timeout=5)  # get remote object
        if response.status_code == 200:  # SUCCESS
            return response.text
        logger.debug(f"GET failed with status {response.status_code} for: {url}")
        return False

    def parse_attributes(self, response_text):  # inherited from Sammyjava
        """parses attributes returned from HTMLParser. Credit to SammyJava"""
        collections = (
            []
        )  # Prevents self.collections collision with self in CollectionsParser
        collection_types = self.collection_types

        class CollectionsParser(HTMLParser):
            """HTMLParser for Collections"""

            def handle_starttag(self, tag, attrs):
                """Feed from HTMLParser"""
                for attr in attrs:  # check each attribute passed in attrs
                    if attr[0] == "href":  # attribute is a URL
                        for collection_type in collection_types:
                            if (
                                f"/{collection_type}/" in attr[1]
                            ):  # collection_type iterator is in URL attr[1]
                                collections.append(
                                    attr[1]
                                )  # add attr[1] to collections for later use in process_collections

        CollectionsParser().feed(response_text)  # populate collections
        self.collections = collections  # set self.collections

    def get_attributes(self, parts):
        """parse parts return url components"""
        gensp = f"{parts[1][:3].lower()}{parts[2][:2].lower()}"  # make gensp
        strain = parts[-2]  # get strain and key information
        return (gensp, strain)

    def process_collections(self, cmds_only, mode):
        """General method to create a jbrowse-components config or populate a blast db using mode"""
        logger = self.logger
        pathlib.Path(self.out_dir).mkdir(parents=True, exist_ok=True)
        for collection_type in self.collection_types:  # for all collections
            for dsfile in self.files.get(
                collection_type, []
            ):  # for all files in all collections
                cmd = ""
                url = self.files[collection_type][dsfile]["url"]
                if not url:  # do not take objects with no defined link
                    continue
                name = self.files[collection_type][dsfile]["name"]
                version = name.split(".")[-1]
                #                dsname = self.files[collection_type][dsfile]["url"].split("/")[-1]
                genus = self.files[collection_type][dsfile]["genus"]
                taxid = self.files[collection_type][dsfile].get("taxid", 0)
                parent = self.files[collection_type][dsfile]["parent"]
                species = self.files[collection_type][dsfile]["species"]
                infraspecies = self.files[collection_type][dsfile]["infraspecies"]
                filetype = url.split(".")[
                    -3
                ]  # get file type from datastore file name filetype.X.gz
                node = {
                    "filename": name,
                    "filetype": filetype,
                    "canonical_type": filetype,
                    "collection_type": collection_type,
                    "url": url,
                    "counts": self.files[collection_type][dsfile].get("counts", None),
                    "busco": self.files[collection_type][dsfile].get("busco", None),
                    "genus": genus,
                    "species": species,
                    "origin": "LIS",
                    "infraspecies": infraspecies,
                    "derived_from": parent,
                }
                # Published metadata captured off the collection README. Without this
                # the node carries no DOI, no taxid and no assembly conventions, and
                # every consumer has to re-fetch the README to get them back.
                for field in NODE_README_FIELDS:
                    value = self.files[collection_type][dsfile].get(field)
                    if value not in (None, "", [], {}, 0):
                        node[field] = value
                self.file_objects.append(node)  # object for DSCensor node
                ### possibly break out next section into methods: blast, jbrowse, then types

                if collection_type == "genomes":  # add genome
                    if mode == "jbrowse":  # for jbrowse
                        cmd = f"jbrowse add-assembly -n {name} --out {os.path.abspath(self.out_dir)}/ -t bgzipFasta --force"
                        cmd += f' --displayName "{genus.capitalize()} {species} {infraspecies} V{version.replace("gnm", "")} {collection_type.capitalize()}" {url}'
                        # print("url: ", url)
                    elif mode == "blast":  # for blast
                        cmd = f"set -o pipefail -o errexit -o nounset; curl {url} | gzip -dc"  # retrieve genome and decompress
                        cmd += f'| makeblastdb -parse_seqids -out {self.out_dir}/{name} -hash_index -dbtype nucl -title "{genus.capitalize()} {species} {infraspecies} V{version.replace("gnm", "")} {collection_type.capitalize()}"'
                        if taxid:
                            cmd += f" -taxid {taxid}"

                if collection_type == "annotations":  # add annotation
                    if mode == "jbrowse":  # for jbrowse
                        if url.endswith(
                            "faa.gz"
                        ):  # only process non faa annotations in jbrowse
                            continue
                        cmd = f"jbrowse add-track -a {parent[0]} --out {os.path.abspath(self.out_dir)}/ --force"
                        cmd += f' -n "{genus.capitalize()} {species} {infraspecies} V{version.replace("ann", "")} {collection_type.capitalize()}" {url}'
                    elif mode == "blast":  # for blast
                        if not url.endswith(
                            "faa.gz"
                        ):  # only process faa annotations in blast
                            continue
                        cmd = f"set -o pipefail -o errexit -o nounset; curl {url} | gzip -dc"  # retrieve genome and decompress
                        cmd += f'| makeblastdb -parse_seqids -out {self.out_dir}/{name} -hash_index -dbtype prot -title "{genus.capitalize()} {species} {infraspecies} V{version.replace("ann", "")} {collection_type.capitalize()}"'
                        if taxid:
                            cmd += f" -taxid {taxid}"

                if collection_type == "genome_alignments":  # add pair-wise paf files
                    if mode == "jbrowse":  # for jbrowse
                        cmd = f"jbrowse add-track --assemblyNames {','.join(parent)} --out {os.path.abspath(self.out_dir)}/ {url} --force"
                        bam_url = self.files[collection_type][dsfile].get(
                            "bam_url", None
                        )
                        if bam_url:
                            bam_name = os.path.basename(bam_url)
                            cmd += f";jbrowse add-track -n {bam_name} --trackId {bam_name} -a {parent[1]}"
                            cmd += f" --out {os.path.abspath(self.out_dir)}/ --indexFile {bam_url}.bai {bam_url} --force"  # add BAM alignment track for genome_alignments
                    elif mode == "blast":  # for blast
                        continue  # Not blastable at the moment

                if collection_type == "expression":  # add bigwigs

                    if mode == "jbrowse":  # for jbrowse
                        if url.endswith("bw"):
                            bw_name = self.files[collection_type][dsfile].get(
                                "name", None
                            )
                            bw_id = bw_name.split(".")[-2:]
                            project_id = ".".join(bw_name.split(".")[1:-2])
                            cmd = f"jbrowse add-track {url} --name {bw_id[0]} --assemblyNames {parent[0]} --category expression,{project_id} --out {os.path.abspath(self.out_dir)} --force"

                    elif mode == "blast":  # for blast
                        continue  # Not blastable at the moment

                # MORE CANONICAL TYPES HERE
                if not cmd:  # continue for null or incomplete objects
                    continue
                if cmds_only:  # output only cmds
                    print(cmd)

                elif subprocess.check_call(
                    cmd, shell=True, executable="/bin/bash"
                ):  # execute cmd and check exit value = 0
                    logger.error(f"Non-zero exit value: {cmd}")

    def populate_jbrowse2(self, out_dir, cmds_only=False):
        """Populate jbrowse2 config object from collected objects"""
        if out_dir:  # set output directory
            self.out_dir = out_dir
        pathlib.Path(self.out_dir).mkdir(parents=True, exist_ok=True)
        self.process_collections(
            cmds_only, "jbrowse"
        )  # process collections for jbrowse-components

    def populate_blast(self, out_dir, cmds_only=False):
        """Populate a BLAST db for genome_main, mrna/mrna_primary and protein/protein_primary"""
        if out_dir:  # set output directory
            self.out_dir = out_dir
        pathlib.Path(self.out_dir).mkdir(parents=True, exist_ok=True)
        self.process_collections(
            cmds_only, "blast"
        )  # process collections for BLAST sequenceserver

    def populate_dscensor(self, out_dir):
        """Populate dscensor nodes for loading into a neo4j database"""
        if out_dir:  # set output directory
            self.out_dir = out_dir
        pathlib.Path(self.out_dir).mkdir(parents=True, exist_ok=True)
        self.process_collections(True, "dscensor")  # process collections for DSCensor
        for node in self.file_objects:  # write all processed objects to node files
            node_out = open(
                f'{self.out_dir}/{node["filename"]}.json', "w", encoding="utf-8"
            )  # file to write node to
            node_out.write(json.dumps(node))
            node_out.close()

    def parse_busco(self, busco_url):
        """Grab BUSCOs from remote busco_url"""
        logger = self.logger
        busco_response = None
        if self.from_github:
            if os.path.isfile(busco_url):
                busco_response = open(busco_url, encoding="utf-8").read()
        else:
            busco_response = self.get_remote(busco_url)
        if not busco_response:
            logger.debug(busco_response)
            return {}
        logger.debug(f"Adding BUSCO: {busco_response}")
        busco_data = json.loads(busco_response)
        logger.debug(busco_data)
        results = busco_data["results"]  # get stats from run for genome
        complete = float(results["Complete"] / 100)
        single_copy = float(results["Single copy"] / 100)
        multi_copy = float(results["Multi copy"] / 100)
        fragmented = float(results["Fragmented"] / 100)
        missing = float(results["Missing"] / 100)
        markers = float(results["n_markers"])
        records = int(results.get("Number of scaffolds", 0))
        contigs = int(results.get("Number of contigs", 0))
        all_bases = int(results.get("Total length", 0))
        n50 = int(results.get("Scaffold N50", 0))
        gap_bases = 0
        if records:
            gap_bases = int(
                all_bases * float(results["Percent gaps"].replace("%", "")) / 100
            )
        else:
            logger.debug(f"No FASTA Stats for: {busco_url}")
        busco_return = {
            "complete_buscos": int(complete * markers),
            "single_copy_buscos": int(single_copy * markers),
            "duplicate_buscos": int(multi_copy * markers),
            "fragmented_buscos": int(fragmented * markers),
            "missing_buscos": int(missing * markers),
            "total_buscos": int(markers),
        }
        genome_return = {
            "contigs": contigs,
            "records": records,
            "N50": n50,
            "allbases": all_bases,
            "gapbases": gap_bases,
            "gaps": contigs - 1,  # this is a hack for now with the gaps value
        }
        gff_return = (
            {}
        )  # add this later if we decide to process gff stats here self.gff3_stats
        if "genomes" in busco_url:
            return {"counts": genome_return, "busco": busco_return}  # return
        if "annotations" in busco_url:
            return {"counts": gff_return, "busco": busco_return}
        return {}

    def add_collections(self, collection_type, genus, species):
        """Adds collection to self.files[collection_type] for later use"""
        logger = self.logger
        logger.debug("in add_collections")
        from_github = self.from_github
        species_url = f"{self.datastore_url}/{genus}/{species}"
        collections_url = f"{species_url}/{collection_type}/"
        collections_dir = None
        collections_response = None
        if from_github:
            species_url = f"{self.from_github}/{genus}/{species}"
            collections_dir = f"{species_url}/{collection_type}/"
            if os.path.isdir(collections_dir):
                collections_response = collections_dir
        else:
            collections_response = self.get_remote(collections_url)
        if not collections_response:  # get remote failed
            logger.debug(collections_response)
            return False
        if collection_type not in self.files:  # add new type
            self.files[collection_type] = {}
        print(
            f"  {collection_type}:", file=self.species_collections_handle
        )  # print collection type in species collections
        self.collections = []  # Set to empty list for use in self.parse_attributes
        if from_github:
            #            for d in os.walk(collections_dir):
            for collection_dir in next(os.walk(collections_dir))[1]:
                collection = "/".join(collections_dir.split("/")[-4:])
                self.collections.append(f"/{collection}{collection_dir}/")
        else:
            self.parse_attributes(
                collections_response
            )  # Feed response from GET to populate collections
        for collection_dir in self.collections:
            parts = collection_dir.split("/")
            #            print(collection_dir, parts)
            #            ['', 'falafel', 'ctc', 'sw', 'LIS-autocontent', 'datastore-metadata', 'Arachis', 'hypogaea', 'genomes', '']
            #            ['', 'Arachis', 'hypogaea', 'genomes', 'BaileyII.gnm1.1JTF', '']
            #            sys.exit(1)
            logger.debug(parts)
            name = parts[4]
            url = ""
            parent = ""
            parts = self.get_attributes(parts)
            lookup = f"{parts[0]}.{'.'.join(name.split('.')[:-1])}"  # reference name in datastructure
            strain_lookup = lookup.split(".")[1]  # the strain for the lookup
            if collection_type == "genomes":  # add parent genome_main files
                ref = ""
                stop = 0
                url = f"{self.datastore_url}{collection_dir}{parts[0]}.{parts[1]}.genome_main.fna.gz"  # genome_main in datastore_url
                fai_url = f"{url}.fai"  # get fai file for jbrowse session construction
                fai_response = self.get_remote(
                    fai_url
                )  # get fai file to build loc from
                if fai_response:  # fai SUCCESS 200
                    ref, stop = fai_response.split("\n")[0].split()[
                        :2
                    ]  # fai field 1\s+2. field 1 is sequence_id field 2 is length
                    logger.debug(f"{ref},{stop}")
                else:  # fai file could not be accessed
                    logger.error(f"No fai file for: {url}")
                    sys.exit(1)

                linear_session = {  # LinearGenomeView object for JBrowse2
                    "views": [
                        {
                            "assembly": lookup,
                            "loc": f"{ref}:1-1000000",  # JBrowse2 does not allow null loc
                            "type": "LinearGenomeView",
                            #                                            "tracks": [
                            #                                                " gff3tabix_genes " ,
                            #                                                " volvox_filtered_vcf " ,
                            #                                                " volvox_microarray " ,
                            #                                                " volvox_cram "
                            #                                            ]
                        }
                    ]
                }
                linear_url = f"{self.jbrowse_url}/?config=config.json&session=spec-{linear_session}"  # build the URL for the resource

                linear_data = {
                    "name": f"JBrowse2 {lookup}",
                    "URL": str(linear_url).replace(
                        "'", "%22"
                    ),  # url encode for .yml file and Jekyll linking
                    "description": "JBrowse2 Linear Genome View",
                }  # the object that will be written into the .yml file

                logger.debug(f"linear data for assembly: {linear_data} \n")

                if strain_lookup not in self.infraspecies_resources:
                    self.infraspecies_resources[strain_lookup] = (
                        []
                    )  # initialize infraspecies list within species
                if self.jbrowse_url:  # dont add data if no jbrowse url set
                    self.infraspecies_resources[strain_lookup].append(linear_data)
                logger.debug(url)
                busco_url = f"{self.datastore_url}{collection_dir}/BUSCO/{parts[0]}.{parts[1]}.busco.fabales_odb10.short_summary.json"
                if from_github:
                    busco_url = f"{self.from_github}/{collection_dir}/BUSCO/{parts[0]}.{parts[1]}.busco.fabales_odb10.short_summary.json"
                logger.debug(busco_url)
                genome_stats = self.parse_busco(busco_url)
                logger.debug(genome_stats)
                if not genome_stats:
                    logger.debug(f"No short summary for: {busco_url}")
                self.files[collection_type][lookup] = {
                    "url": url,
                    "name": lookup,
                    "parent": [parent],
                    "genus": genus,
                    "species": species,
                    "infraspecies": strain_lookup,
                    "taxid": 0,
                    "busco": genome_stats.get("busco"),
                    "counts": genome_stats.get("counts"),
                }  # add type and lookup for object with labels, stats and buscos
                logger.debug(self.files[collection_type][lookup])
            ###
            elif (
                collection_type == "annotations"
            ):  # add gff3 annotation files and protein/protein_primary. genome_main parent
                genome_lookup = ".".join(lookup.split(".")[:-1])  # genome parent prefix
                #                self.files["genomes"][genome_lookup]["url"]
                parent = genome_lookup
                url = f"{self.datastore_url}{collection_dir}{parts[0]}.{parts[1]}.gene_models_main.gff3.gz"
                busco_url = f"{self.datastore_url}{collection_dir}/BUSCO/{parts[0]}.{parts[1]}.busco.fabales_odb10.short_summary.json"
                if from_github:
                    busco_url = f"{self.from_github}/{collection_dir}/BUSCO/{parts[0]}.{parts[1]}.busco.fabales_odb10.short_summary.json"
                logger.debug(busco_url)
                annotation_stats = self.parse_busco(busco_url)
                logger.debug(annotation_stats)
                if not annotation_stats:
                    logger.debug(f"No short summary for: {busco_url}")
                self.files[collection_type][lookup] = {  # gene_models_main
                    "url": url,
                    "name": lookup,
                    "parent": [parent],
                    "genus": genus,
                    "species": species,
                    "infraspecies": strain_lookup,
                    "taxid": 0,
                    "busco": annotation_stats.get("busco"),
                    "counts": annotation_stats.get("counts"),
                }  # add type and url
                logger.debug(self.files[collection_type][lookup])
                protprimary_url = f"{self.datastore_url}{collection_dir}{parts[0]}.{parts[1]}.protein_primary.faa.gz"
                protprimary_response = self.head_remote(protprimary_url)
                if protprimary_response:
                    protprimary_lookup = f"{lookup}.protein_primary"
                    self.files[collection_type][protprimary_lookup] = (
                        {  # protein_primary
                            "url": protprimary_url,
                            "name": protprimary_lookup,
                            "parent": [parent],
                            "genus": genus,
                            "species": species,
                            "infraspecies": strain_lookup,
                            "taxid": 0,
                        }
                    )
                else:
                    logger.debug(
                        f"protein_primary failed:{protprimary_url}, {protprimary_response}"
                    )

                protein_url = f"{self.datastore_url}{collection_dir}{parts[0]}.{parts[1]}.protein.faa.gz"
                protein_response = self.head_remote(protein_url)
                if protein_response:
                    protein_lookup = f"{lookup}.protein"
                    self.files[collection_type][protein_lookup] = {  # all proteins
                        "url": protein_url,
                        "name": protein_lookup,
                        "parent": [parent],
                        "genus": genus,
                        "species": species,
                        "infraspecies": strain_lookup,
                        "taxid": 0,
                    }
                else:
                    logger.debug(f"protein failed:{protein_url}, {protein_response}")
            ###
            #            elif collection_type == "synteny":  # DEPRICATED?
            #                checksum_url = f"{self.datastore_url}{collection_dir}CHECKSUM.{parts[1]}.md5"
            #                checksum_response = requests.get(checksum_url)
            #                if checksum_response.status_code == 200:
            #                    continue
            #                else:  # CheckSum FAILURE
            #                    logger.debug(
            #                        f"GET Failed for checksum {checksum_response.status_code} {checksum_url}"
            #                    )
            ###
            elif (
                collection_type == "genome_alignments"
            ):  # Synteny after the new changes. Parent is a tuple with both genome_main files
                checksum_url = (
                    f"{self.datastore_url}{collection_dir}CHECKSUM.{parts[1]}.md5"
                )
                checksum_response = None
                if from_github:
                    checksum_response = open(
                        f"{self.from_github}/{collection_dir}CHECKSUM.{parts[1]}.md5",
                        encoding="utf-8",
                    ).read()
                else:
                    checksum_response = self.get_remote(checksum_url)
                logger.debug(checksum_response)
                if checksum_response:  # checksum SUCCESS 200
                    for line in checksum_response.split("\n"):
                        logger.debug(line)
                        fields = line.split()
                        if fields:  # process if fields exists
                            if fields[1].endswith("paf.gz"):  # get paf file
                                paf_lookup = fields[1].replace(
                                    "./", ""
                                )  # get paf file to load will start with ./
                                logger.debug(paf_lookup)
                                paf_url = f"{self.datastore_url}{collection_dir}{paf_lookup}"  # where the paf file is in the datastore
                                paf_parts = paf_lookup.split(
                                    "."
                                )  # split the paf file name into parts delimited by '.'
                                parent1 = ".".join(
                                    paf_parts[:3]
                                )  # parent 1 in pair-wise alignment
                                parent2 = ".".join(
                                    paf_parts[4:7]
                                )  # parent 2 in pair-wise alignment
                                self.files[collection_type][paf_lookup] = {
                                    "url": paf_url,
                                    "name": paf_lookup,
                                    "parent": [parent2, parent1],
                                    "genus": genus,
                                    "species": species,
                                    "infraspecies": strain_lookup,
                                    "taxid": 0,
                                    "bam_url": paf_url.replace("paf.gz", "bam"),
                                }
                                logger.debug(self.files[collection_type][paf_lookup])
                                dotplot_view = {  # session object for jbrowse2 dotplot view populate below with parent1 and parent2
                                    "views": [
                                        {
                                            "type": "DotplotView",
                                            "views": [
                                                {"assembly": parent1},
                                                {"assembly": parent2},
                                            ],
                                            "tracks": [paf_lookup.replace(".gz", "")],
                                        }
                                    ]
                                }
                                dotplot_url = f"{self.jbrowse_url}/?config=config.json&session=spec-{dotplot_view}"  # build the URL for the resource
                                dotplot_data = {
                                    "name": f"JBrowse2 {paf_lookup}",
                                    "URL": str(dotplot_url).replace(
                                        "'", "%22"
                                    ),  # url encode for .yml file and Jekyll linking
                                    "description": "JBrowse2 Dotplot View",
                                }  # the object that will be written into the .yml file
                                if strain_lookup not in self.infraspecies_resources:
                                    self.infraspecies_resources[strain_lookup] = (
                                        []
                                    )  # initialize infraspecies list within species
                                if (
                                    self.jbrowse_url
                                ):  # dont add data if no jbrowse url set
                                    self.infraspecies_resources[strain_lookup].append(
                                        dotplot_data
                                    )  # add data for later writing in resources
            ###
            elif collection_type == "expression":  # add parent expr files
                ref = ""
                # Synteny after the new changes. Parent is a tuple with both genome_main files
                checksum_url = (
                    f"{self.datastore_url}{collection_dir}CHECKSUM.{parts[1]}.md5"
                )
                # print('\nchecksum_url: ', checksum_url,'\n')
                checksum_response = None
                if from_github:
                    checksum_response = open(
                        f"{self.from_github}/{collection_dir}CHECKSUM.{parts[1]}.md5",
                        encoding="utf-8",
                    ).read()
                else:
                    checksum_response = self.get_remote(checksum_url)
                logger.debug(checksum_response)
                # print('\nchecksum_response: ', checksum_response,'\n')
                if checksum_response:  # checksum SUCCESS 200
                    for line in checksum_response.split("\n"):
                        logger.debug(line)
                        fields = line.split()
                        if fields:  # process if fields exists
                            if fields[1].endswith("bw"):  # get bw file
                                ref = ""
                                stop = 0
                                bw_lookup = fields[1].replace(
                                    "./", ""
                                )  # get bw file to load will start with ./
                                logger.debug(bw_lookup)
                                bw_url = f"{self.datastore_url}{collection_dir}{bw_lookup}"  # where the bigwig file is in the datastore
                                bw_parts = bw_lookup.split(
                                    "."
                                )  # split the bw file name into parts delimited by '.'
                                genome_lookup = ".".join(
                                    lookup.split(".")[:-3]
                                )  # genome parent prefix
                                #               self.files["genomes"][genome_lookup]["url"]
                                parent = genome_lookup
                                self.files[collection_type][bw_lookup] = {
                                    "url": bw_url,
                                    "name": bw_lookup,
                                    "genus": genus,
                                    "parent": [parent],
                                    "species": species,
                                    "infraspecies": strain_lookup,
                                    "taxid": 0,
                                }
                                logger.debug(self.files[collection_type][bw_lookup])

                                # url =  f"{self.datastore_url}{collection_dir}{parts[0]}.{parts[1]}.genome_main.fna.gz"  # genome_main in datastore_url
                                url = self.files["genomes"][parent]["url"]
                                fai_url = f"{url}.fai"  # get fai file for jbrowse session construction
                                fai_response = self.get_remote(
                                    fai_url
                                )  # get fai file to build loc from
                                if fai_response:  # fai SUCCESS 200
                                    ref, stop = fai_response.split("\n")[0].split()[
                                        :2
                                    ]  # fai field 1\s+2. field 1 is sequence_id field 2 is length
                                    logger.debug(f"{ref},{stop}")
                                else:  # fai file could not be accessed
                                    logger.error(f"No fai file for: {url}")
                                    sys.exit(1)

                                linear_session = {  # LinearGenomeView object for JBrowse2
                                    "views": [
                                        {
                                            "assembly": parent,
                                            # sequence is currently hardcoded, don't know how "ref" works for genomes
                                            "loc": f"{ref}:1-{stop}",  # JBrowse2 does not allow null loc
                                            "type": "LinearGenomeView",
                                            "tracks": [
                                                ".".join(bw_lookup.split(".")[:-1])
                                            ],
                                            # ["glyma.Wm82.gnm6.ann1.expr.mixed.Kour_Boone_2014.Clark_defective"]
                                            #                                                " gff3tabix_genes " ,
                                            #                                                " volvox_filtered_vcf " ,
                                            #                                                " volvox_microarray " ,
                                            #                                                " volvox_cram "
                                            #                                            ]
                                        }
                                    ]
                                }
                                linear_url = f"{self.jbrowse_url}/?config=config.json&session=spec-{linear_session}"  # build the URL for the resource

                                linear_data = {
                                    "name": f"JBrowse2 {parent}",
                                    "URL": str(linear_url).replace(
                                        "'", "%22"
                                    ),  # url encode for .yml file and Jekyll linking
                                    "description": "JBrowse2 Linear Genome View",
                                }  # the object that will be written into the .yml file

                                if strain_lookup not in self.infraspecies_resources:
                                    self.infraspecies_resources[strain_lookup] = (
                                        []
                                    )  # initialize infraspecies list within species
                                if (
                                    self.jbrowse_url
                                ):  # dont add data if no jbrowse url set
                                    self.infraspecies_resources[strain_lookup].append(
                                        linear_data
                                    )
                                logger.debug(f"linear data for bw: {linear_data} \n")
            ###

            readme_url = f"{self.datastore_url}/{collection_dir}README.{name}.yml"  # species collection readme
            readme_response = None
            if from_github:
                github_readme = f"{self.from_github}/{collection_dir}README.{name}.yml"
                if os.path.isfile(github_readme):
                    readme_response = open(github_readme, encoding="utf-8").read()
            else:
                readme_response = self.get_remote(readme_url)
            if readme_response:  # readme get success
                readme = yaml.load(readme_response, Loader=yaml.FullLoader)
                logger.debug(readme)
                synopsis = readme["synopsis"]
                # Carry the published metadata onto every file entry from this
                # collection, not just the collection-level one: a collection yields
                # several entries ("<lookup>", "<lookup>.protein", ...) and each becomes
                # its own node. Anything not copied here is silently lost -- which is
                # how publication_doi and taxid used to disappear.
                metadata = {}
                for field in NODE_README_FIELDS:
                    value = readme.get(field)
                    if value not in (None, "", [], {}):
                        metadata[field] = value
                matched = 0
                for key in self.files[collection_type]:
                    if key == lookup or key.startswith(f"{lookup}."):
                        self.files[collection_type][key].update(metadata)
                        matched += 1
                if not matched:
                    logger.debug(f"{lookup} not in {self.files[collection_type]}")
                print(
                    f"    - collection: {name}",
                    file=self.species_collections_handle,
                )
                print(
                    f'      synopsis: "{synopsis}"',
                    file=self.species_collections_handle,
                )
            else:  # get failed for
                logger.debug(f"GET Failed for README {readme_url}")
        return True

    def process_species(self, genus, species):
        """Process species and genus from genus_description object"""
        logger = self.logger
        logger.debug("in process_species")
        from_github = self.from_github
        if from_github:
            logger.info(f"Searching {self.from_github} for: {genus} {species}")
        else:
            logger.info(f"Searching {self.datastore_url} for: {genus} {species}")
        species_url = f"{self.datastore_url}/{genus}/{species}"
        self.infraspecies_resources = {}
        print(f"- name: {species}", file=self.species_collections_handle)

        for (
            collection_type
        ) in (
            self.collection_types
        ):  # iterate through collections found in the datastore
            self.add_collections(collection_type, genus, species)

        species_description_url = f"{species_url}/about_this_collection/description_{genus}_{species}.yml"  # parse for strain resources
        logger.debug(species_description_url)  # get species description url
        species_description_response = None
        if from_github:
            species_description_response = open(
                f"{from_github}/{genus}/{species}/about_this_collection/description_{genus}_{species}.yml",
                encoding="utf-8",
            ).read()
        else:
            species_description_response = self.get_remote(species_description_url)
        if (
            species_description_response
        ):  # Read species description yml and add jbrowse resources to "strains"
            species_description = yaml.load(
                species_description_response, Loader=yaml.FullLoader
            )  # load the yaml from the datastore for species
            count = 0
            for strain in species_description[
                "strains"
            ]:  # iterate through all strains in this species description
                if (
                    strain["identifier"] in self.infraspecies_resources
                ):  # add to this strain
                    if species_description["strains"][count].get(
                        "resources", None
                    ):  # this strain has resources
                        for resource in self.infraspecies_resources[
                            strain["identifier"]
                        ]:  # append all the resources to the existing
                            species_description["strains"][count]["resources"].append(
                                resource
                            )
                    else:
                        species_description["strains"][count]["resources"] = (
                            self.infraspecies_resources[strain["identifier"]]
                        )  # set resources
                count += 1  # keep track of how many "strains" we have seen
            self.species_descriptions.append(species_description)

    def process_taxon(self, taxon):
        """Retrieve and output collections for jekyll site"""
        logger = self.logger
        from_github = self.from_github
        logger.debug(f"in process_taxon {from_github}")
        if not "genus" in taxon:  # genus required for all taxon
            logger.error(f"Genus not found for: {taxon}")
            sys.exit(1)
        genus = taxon["genus"]
        genus_description_url = f"{self.datastore_url}/{genus}/GENUS/about_this_collection/description_{genus}.yml"  # genus description to be read
        genus_description_response = None
        if from_github:  # if build locally from github clone of datastore-metadata
            genus_description_response = open(
                f"{self.from_github}/{genus}/GENUS/about_this_collection/description_{genus}.yml",
                encoding="utf-8",
            ).read()
        else:
            genus_description_response = self.get_remote(genus_description_url)
        logger.debug(genus_description_response)
        self.genus_resources_handle = None  # yaml file to write for genus resources
        self.species_resources_handle = None  # yaml file to write for species resources
        self.species_collections_handle = (
            None  # yaml file to write for species collections
        )
        if genus_description_response:  # Genus Description yml 200 SUCCESS
            species_collections_filename = None
            self.species_descriptions = []  # null for current taxon genus
            genus_description = yaml.load(
                genus_description_response, Loader=yaml.FullLoader
            )  # load yml into python object
            collection_dir = f"{os.path.abspath(self.out_dir)}/{genus}"
            pathlib.Path(collection_dir).mkdir(
                parents=True, exist_ok=True
            )  # make output dirs if they dont exist
            genus_resources_filename = f"{collection_dir}/genus_resources.yml"  # local file to write genus resources
            species_resources_filename = f"{collection_dir}/species_resources.yml"  # local file to write species resources
            species_collections_filename = f"{collection_dir}/species_collections.yml"  # local file to write collections
            self.species_collections_handle = open(
                species_collections_filename, "w", encoding="utf-8"
            )
            self.genus_resources_handle = open(
                genus_resources_filename, "w", encoding="utf-8"
            )
            self.species_resources_handle = open(
                species_resources_filename, "w", encoding="utf-8"
            )
            collection_string = "---\nspecies:"
            print("---", file=self.genus_resources_handle)  # write genus resources
            yaml.dump(
                genus_description, self.genus_resources_handle
            )  # dump full description of genus
            print(
                collection_string, file=self.species_collections_handle
            )  # write species collection
            print(
                collection_string, file=self.species_resources_handle
            )  # write species resources

            for species in genus_description[
                "species"
            ]:  # iterate through all species in the genus
                self.process_species(
                    genus, species
                )  # process species and genus to populate self.species_descriptions

            yaml.dump(
                self.species_descriptions, self.species_resources_handle
            )  # dump species_resources.yml locally with all self.species_descriptions from genus
        if self.genus_resources_handle:  # close genus resources
            self.genus_resources_handle.close()
        if self.species_resources_handle:  # close species resources
            self.species_resources_handle.close()
        if self.species_collections_handle:  # close species collections
            self.species_collections_handle.close()

    def discover_genera(self):
        """Return every genus in the datastore clone, sorted.

        A directory counts as a genus only when it holds a GENUS description
        file, which excludes non-taxonomic tops such as LEGUMES and .github.
        """
        root = pathlib.Path(self.from_github)
        candidates = sorted(
            entry.name
            for entry in root.iterdir()
            if entry.is_dir() and not entry.name.startswith(".")
        )
        genera = [
            genus
            for genus in candidates
            if (
                root
                / genus
                / "GENUS"
                / "about_this_collection"
                / f"description_{genus}.yml"
            ).is_file()
        ]
        skipped = sorted(set(candidates) - set(genera))
        if skipped:
            self.logger.debug(f"Not genus directories, skipped: {skipped}")
        self.logger.info(f"Discovered {len(genera)} genera: {genera}")
        return genera

    def load_taxa(self, target=None):
        """Return the taxa to process, from a taxon list file or by discovery."""
        if target and os.path.isdir(target):
            self.logger.warning(
                f"Taxon list {target} is a directory, discovering genera instead"
            )
            target = None
        if target:
            self.logger.info(f"Reading taxon list: {target}")
            with open(target, "r", encoding="utf-8") as handle:
                return yaml.load(handle.read(), Loader=yaml.FullLoader)
        if not self.from_github:
            self.logger.error(
                "A taxon list is required when not reading a datastore-metadata clone"
            )
            sys.exit(1)
        self.logger.info("No taxon list given, discovering all genera")
        return [{"genus": genus} for genus in self.discover_genera()]

    def parse_collections(
        self, from_github="./datastore-metadata", target=None
    ):  # refactored from SammyJava
        """Retrieve and output collections for jekyll site.

        from_github is the datastore-metadata clone to read. With no target,
        every genus found there is processed; pass a taxon list yml to restrict
        the run to a subset.
        """
        if from_github:  # set to None if empty dir
            self.from_github = os.path.abspath(from_github)
        self.logger.debug(f"THIS IS GITHUB: {self.from_github}")
        for taxon in self.load_taxa(target):
            self.process_taxon(taxon)  # process taxon object


if __name__ == "__main__":
    parser = ProcessCollections()
    parser.parse_collections()
