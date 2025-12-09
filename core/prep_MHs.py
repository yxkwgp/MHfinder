#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microhaplotype preparation functionality module
"""

import os
from random import choice
from .utils import run_plink_command


class MicrohaplotypeGenerator:
    """Microhaplotype preparation class"""
    
    def _is_not_subset_any(self, checked_set, sets_list):
        """Check if a set is not a subset of any set in the list"""
        for element in sets_list:
            if checked_set.issubset(element):
                return False
        return True
    
    def _find_cliques(self, unmatch):
        """Bron Kerbosch method to find maximal cliques"""
        cliques = []
        def bron_kerbosch(R=set(), P=set(), X=set()):
            if not P and not X:
                if len(R) > 1:
                    cliques.append(R)
                return
            u = choice(list(P.union(X)))
            for v in P.difference(unmatch[u]):
                newR = R.copy()
                newR.add(v)
                newP = P.intersection(unmatch[v])
                newX = X.intersection(unmatch[v])
                bron_kerbosch(newR, newP, newX)
                P.remove(v)
                X.add(v)
            return
        bron_kerbosch(P=set(range(len(unmatch))))
        return cliques
    
    def prepare_microhaplotypes(self, input_prefix, snplist, length_threshold, 
                               min_site_number, r2_threshold, output):
        """
        Generate microhaplotypes
        
        Parameters:
            input_prefix (str): PLINK binary file prefix
            snplist (str): Filtered snplist file
            length_threshold (int): Length threshold (distance threshold between two endpoint SNPs of microhaplotype)
            min_site_number (int): Minimum number of SNPs in microhaplotype
            r2_threshold (float): Strong linkage threshold (r2 between any two SNPs in microhaplotype must be below this value)
            output (str): Final output filename
            
        Returns:
            str: Generated microhaplotype file path
        """
        # Parameter validation
        if not os.path.exists(f"{input_prefix}.bed"):
            raise FileNotFoundError(f"PLINK file not found: {input_prefix}.bed")
        if not os.path.exists(snplist):
            raise FileNotFoundError(f"SNP list file not found: {snplist}")
        if length_threshold <= 0:
            raise ValueError("length_threshold must be greater than 0")
        if min_site_number < 2:
            raise ValueError("min_site_number must be greater than or equal to 2")
        if r2_threshold < 0 or r2_threshold > 1:
            raise ValueError("r2_threshold must be between 0 and 1")
        
        print(f"Starting microhaplotype generation, length threshold: {length_threshold}, minimum sites: {min_site_number}")
        
        # plink calculate linkage disequilibrium
        length_threshold_kb = length_threshold / 1000
        cmd = f'plink --bfile {input_prefix} --extract {snplist} --r2 --ld-window 10000 --ld-window-kb {length_threshold_kb} --ld-window-r2 0 --out {input_prefix}'
        run_plink_command(cmd)

        ld_file = f'{input_prefix}.ld'
        if not os.path.exists(ld_file):
            raise RuntimeError(f"Linkage disequilibrium file not generated: {ld_file}")
        
        outfile = open(output, 'w')
        outfile.write('chr|SNP1-SNP2...SNPn\n')
        infile1 = open(ld_file, 'r')
        infile2 = open(ld_file, 'r')
        infile1.readline()
        infile2.readline()
        
        # First line of ld file
        anchor_line = infile1.readline().strip().split()
        site_group = []
        if float(anchor_line[6]) < r2_threshold and int(anchor_line[1]) + length_threshold > int(anchor_line[4]):
            site_group.append(anchor_line)
        
        # Initialize ld records
        ld_rem_pos = []
        ld_rem_pair = []
        ld_rem_r2 = []
        ld_line = infile2.readline().strip().split()
        
        # Initialize MH prototype library
        mhs_pri_rem = []
        
        # Process line by line
        for line in infile1:
            line = line.strip().split()
            if line[2] == anchor_line[2]:
                if float(line[6]) < r2_threshold and int(line[1]) + length_threshold > int(line[4]):
                    site_group.append(line)
            else:
                if len(site_group) >= min_site_number - 1:
                    mh_pri = set([site_group[0][2]] + [site_group[i][5] for i in range(len(site_group))])
                    if self._is_not_subset_any(mh_pri, mhs_pri_rem):
                        mhs_pri_rem.append(mh_pri)
                        if len(mhs_pri_rem) > length_threshold:
                            mhs_pri_rem = mhs_pri_rem[1:]
                        
                        front_site_info = [site_group[0][i] for i in [2, 0, 1]]  # [rsID, chr, bp]
                        headless_mh_waiting_check = [site_group[i][5] for i in range(len(site_group))]
                        
                        while ld_line and len(ld_line) > 0 and ld_line[0] != front_site_info[1]:
                            next_line = infile2.readline()
                            if not next_line:  # End of file
                                break
                            ld_line = next_line.strip().split()
                        
                        while ld_line and len(ld_line) > 0 and ld_line[0] == front_site_info[1] and int(ld_line[1]) < int(front_site_info[2]) + length_threshold:
                            ld_rem_pos.append(int(ld_line[1]))
                            ld_rem_pair.append(f'{ld_line[2]}-{ld_line[5]}')
                            ld_rem_r2.append(float(ld_line[6]))
                            next_line = infile2.readline()
                            if not next_line:  # End of file
                                break
                            ld_line = next_line.strip().split()
                        
                        pair_index = 0
                        while pair_index < len(ld_rem_pos) and ld_rem_pos[pair_index] < int(front_site_info[2]):
                            pair_index += 1
                        
                        ld_rem_pos = ld_rem_pos[pair_index:]
                        ld_rem_pair = ld_rem_pair[pair_index:]
                        ld_rem_r2 = ld_rem_r2[pair_index:]
                        
                        unmatch = [set() for i in headless_mh_waiting_check]
                        for m in range(len(headless_mh_waiting_check) - 1):
                            for n in range(m + 1, len(headless_mh_waiting_check)):
                                pair = f'{headless_mh_waiting_check[m]}-{headless_mh_waiting_check[n]}'
                                if pair in ld_rem_pair:
                                    r2 = ld_rem_r2[ld_rem_pair.index(pair)]
                                    if r2 < r2_threshold:
                                        unmatch[m].add(n)
                                        unmatch[n].add(m)
                        
                        cliques = self._find_cliques(unmatch)
                        for clique in cliques:
                            if len(clique) >= min_site_number - 1:
                                clique = list(clique)
                                clique.sort()
                                headless_mh_formal = [headless_mh_waiting_check[i] for i in clique]
                                outfile.write(f'{front_site_info[1]}|{front_site_info[0]}-{"-".join(headless_mh_formal)}\n')
                
                # Update anchorLine
                site_group = []
                anchor_line = line[:]
                if float(anchor_line[6]) < r2_threshold:
                    site_group.append(anchor_line)
        
        infile1.close()
        infile2.close()
        outfile.close()
        
        print(f"Microhaplotype generation completed, output file: {output}")
        return output 