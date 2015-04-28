from collections import defaultdict
import sys

import numpy as np
import gffutils
import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def miso_exon_to_gencode_exon(exon):
    """Convert a single miso exon to gencode gffutils database exon id

    >>> # Skipped exon (SE) or Mutually exclusive exon (MXE) ID
    >>> miso_exon_to_gencode_exon('chr2:9624561:9624679:+')
    'exon:chr2:9624561-9624679:+'
    >>> # Alt 5' splice site (pick first of alternative)
    >>> miso_exon_to_gencode_exon('chr15:42565276:42565087|42565161:-')
    'exon:chr15:42565276-42565087:-'
    >>> # Alt 3' splice site (pick first of alternative)
    >>> miso_exon_to_gencode_exon('chr2:130914199|130914248:130914158:-')
    'exon:chr2:130914199-130914158:-'
    >>> # Retained intron: start, stop separated by '-' instead of ':'
    >>> miso_exon_to_gencode_exon('chr1:906259-906386:+')
    'exon:chr1:906259-906386:+'
    """
    return 'exon:{}:{}-{}:{}'.format(*miso_exon_to_coords(exon))


def miso_id_to_exon_ids(miso_id):
    """Convert a MISO-style alternative event ID to a gffutils exon id of
    all exons in all possible transcripts

    Split on the pipe ("|") to account for Alt 5'/3' splice site events

    # A skipped exon (SE) ID
    >>> miso_id_to_exon_ids('chr2:9624561:9624679:+@chr2:9627585:9627676:+@chr2:9628276:9628591:+')
    ['exon:chr2:9624561-9624679:+', 'exon:chr2:9627585-9627676:+', 'exon:chr2:9628276-9628591:+']
    >>> # A mutually exclusive (MXE) ID
    >>> miso_id_to_exon_ids('chr16:89288500:89288591:+@chr16:89289565:89289691:+@chr16:89291127:89291210:+@chr16:89291963:89292039:+')
    ['exon:chr16:89288500-89288591:+', 'exon:chr16:89289565-89289691:+', 'exon:chr16:89291127-89291210:+', 'exon:chr16:89291963-89292039:+']
    >>> # An Alt 5' splice site (A5SS) ID
    >>> miso_id_to_exon_ids("chr15:42565276:42565087|42565161:-@chr15:42564261:42564321:-")
    ['exon:chr15:42565276-42565087:-', 'exon:chr15:42565276-42565161:-', 'exon:chr15:42564261-42564321:-']
    >>> # An Alt 3' splice site (A3SS) ID
    >>> miso_id_to_exon_ids('chr2:130914824:130914969:-@chr2:130914199|130914248:130914158:-')
    ['exon:chr2:130914824-130914969:-', 'exon:chr2:130914199-130914158:-', 'exon:chr2:130914248-130914158:-']
    >>> # A retained intron (RI) ID
    >>> miso_id_to_exon_ids('chr1:906066-906138:+@chr1:906259-906386:+')
    ['exon:chr1:906066-906138:+', 'exon:chr1:906259-906386:+', 'exon:chr1:906066-906386:+']
    """
    return map(miso_exon_to_gencode_exon, miso_id.split('@'))


def miso_exon_to_coords(exon):
    """Convert a miso exon to gffutils coordinates

    >>> miso_exon_to_coords('chr2:130914824:130914969:-')
    ('chr2', '130914824', '130914969', '-')
    >>> # Alt 5' SS - pick the first of the alternative ends
    >>> miso_exon_to_coords('chr15:42565276:42565087|42565161:-')
    ('chr15', '42565276', '42565087', '-')
    >>> # Alt 3' SS - pick the first of the alternative starts
    >>> miso_exon_to_coords('chr2:130914199|130914248:130914158:-')
    ('chr2', '130914199', '130914158', '-')
    >>> # Retained intron
    >>> miso_exon_to_coords('chr1:906066-906138:+')
    ('chr1', '906066', '906138', '+')
    """
    strand = exon[-1]
    coords = map(lambda x: x.split('|')[0],
                 exon.split(':'))
    if '-' in coords[1]:
        start, stop = coords[1].split('-')
        coords = coords[0], start, stop, strand
    return coords[0], coords[1], coords[2], strand


def convert_miso_ids_to_everything(miso_ids, db,
                                   event_type,
                                   out_dir):
    """Given a list of miso IDs and a gffutils database, pull out the
    ensembl/gencode/gene name/gene type/transcript names, and write files
    into the out directory. Does not return a value.

    Parameters
    ----------
    miso_ids : list of str
        Miso ids to convert
    db : gffutils.FeatureDB
        gffutils feature database created from a gtf file
    event_type : str
        The type of splicing event. This is used for naming only
    out_dir : str
        Where to write the files to.
    """
    out_dir = out_dir.rstrip('/')
    event_type = event_type.lower()

    miso_to_ensembl = {}
    miso_to_gencode = {}
    miso_to_gene_name = {}
    miso_to_gene_type = {}
    miso_to_ensembl_transcript = {}
    miso_to_gencode_transcript = {}

    ensembl_to_miso = defaultdict(list)
    gencode_to_miso = defaultdict(list)
    gene_name_to_miso = defaultdict(list)

    n_miso_ids = len(miso_ids)
    sys.stdout.write('Converting {} {} miso ids using {} gffutils database '
                     'into {}.'.format(n_miso_ids, event_type, str(db),
                                       out_dir))

    for i, miso_id in enumerate(miso_ids):
        if i % 100 == 0:
            sys.stdout.write('On {}/{} {} miso ids'.format(i, n_miso_ids,
                                                           event_type))

        exons = miso_id_to_exon_ids(miso_id)

        gencode = set([])
        ensembl = set([])
        gene_name = set([])
        gene_type = set([])
        gencode_transcript = set([])
        ensembl_transcript = set([])
        for e in exons:
            try:
                exon = db[e]
                gencode.update(exon.attributes['gene_id'])
                ensembl.update(
                    map(lambda x: x.split('.')[0],
                        exon.attributes['gene_id']))
                gene_name.update(exon.attributes['gene_name'])
                gene_type.update(exon.attributes['gene_type'])
                gencode_transcript.update(exon.attributes['transcript_id'])
                ensembl_transcript.update(
                    map(lambda x: x.split('.')[0], exon.attributes[
                        'transcript_id']))
            except gffutils.FeatureNotFoundError:
                try:
                    #  not an exon, look for any overlapping transcripts here
                    prefix, chrom, startstop, strand = e.split(':')
                    start, stop = startstop.split('-')
                    transcripts = list(db.features_of_type('transcript',
                                                           strand=strand,
                                                           limit=(
                                                           chrom, int(start),
                                                           int(stop))))
                    # if there are overlapping transcripts...
                    if len(transcripts) == 0:
                        continue
                    else:
                        for transcript in transcripts:
                            gencode.update(transcript.attributes['gene_id'])
                            ensembl.update(
                                map(lambda x: x.split('.')[0],
                                    transcript.attributes['gene_id']))
                            gene_name.update(transcript.attributes['gene_name'])
                            gene_type.update(transcript.attributes['gene_type'])
                            gencode_transcript.update(
                                transcript.attributes['transcript_id'])
                            ensembl_transcript.update(
                                map(lambda x: x.split('.')[0],
                                    transcript.attributes['transcript_id']))
                except:
                    continue
        if len(gencode) > 0:

            for ens in ensembl:
                ensembl_to_miso[ens].append(miso_id)
            for g in gene_name:
                gene_name_to_miso[g].append(miso_id)
            for g in gencode:
                gencode_to_miso[g].append(miso_id)

            ensembl = ','.join(ensembl)

            gencode = ','.join(gencode)
            gene_name = ','.join(gene_name)
            gene_type = ','.join(gene_type)

            gencode_transcript = ','.join(gencode_transcript)
            ensembl_transcript = ','.join(ensembl_transcript)

            miso_to_gencode[miso_id] = gencode
            miso_to_ensembl[miso_id] = ensembl
            miso_to_gene_name[miso_id] = gene_name
            miso_to_gene_type[miso_id] = gene_type

            miso_to_gencode_transcript[miso_id] = gencode_transcript
            miso_to_ensembl_transcript[miso_id] = ensembl_transcript

        else:
            miso_to_gencode[miso_id] = np.nan
            miso_to_ensembl[miso_id] = np.nan
            miso_to_gene_name[miso_id] = np.nan
            miso_to_gene_type[miso_id] = np.nan

    miso_tos = {'gencode_gene': miso_to_gencode,
                'ensembl_gene': miso_to_ensembl,
                'gene_name': miso_to_gene_name,
                'gene_type': miso_to_gene_type,
                'gencode_transcript': miso_to_gencode_transcript,
                'ensembl_transcript': miso_to_ensembl_transcript}

    to_misos = {'ensembl_gene': ensembl_to_miso,
                'gene_name': gene_name_to_miso,
                'gencode_gene': gencode_to_miso}

    for name, d in miso_tos.iteritems():
        df = pd.DataFrame.from_dict(d, orient='index')
        df.index.name = 'event_name'
        df.columns = [name]
        tsv = '{}/miso_{}_to_{}.tsv'.format(out_dir, event_type, name)
        df.to_csv(tsv, sep='\t')
        sys.stdout.write('Wrote {}\n'.format(tsv))

    for name, d in to_misos.iteritems():
        tsv = '{}/{}_to_miso_{}.tsv'.format(out_dir, name, event_type)
        with open(tsv, 'w') as f:
            for k, v in d.iteritems():
                f.write('{}\t{}\n'.format(k, '\t'.join(v)))
        sys.stdout.write('Wrote {}\n'.format(tsv))

def splice_type_exons(splice_type, exons):
    """Get exons corresponding to a particular isoform of a splice type

    For SE:
        isoform1: exon1, exon2
        isoform2: exon1, exon2, exon3
    for MXE:
        isoform1: exon1, exon3, exon4
        isoform2: exon1, exon2, exon4

    Parameters
    ----------
    splice_type : 'SE' | 'MXE'
        String specifying the splice type. Currently only SE (skipped exon) and
        MXE (mutually exclusive exon) are supported
    exons : list
        List of exons or CDS's (ids, strings, you name it) in the exact order
        of the splice type, e.g. (exon1_id, exon2_id, exon3_id) for SE

    Returns
    -------
    isoform1_exons : tuple
        Tuple of exons corresponding to isoform 1
    isoform2_exons : tuple
        Tuple of exons corresponding to isoform 2
    """
    isoform1, isoform2 = None, None
    if splice_type == 'SE':
        isoform1 = exons[0], exons[2]
        isoform2 = exons[0], exons[1], exons[2]
    elif splice_type == 'MXE':
        isoform1 = exons[0], exons[1], exons[3]
        isoform2 = exons[0], exons[2], exons[3]
    return isoform1, isoform2
    # if isoform == 1:
    #     return isoform1
    # elif isoform == 2:
    #     return isoform2


def splice_type_isoforms(splice_type, transcripts):
    """Get transcripts corresponding to isoform1 or isoform2 of a splice type

    Parameters
    ----------

    transcripts : list
        List of gffutils iterators of transcripts, where the position indicates
        the exon for which this transcript corresponds.
        For example, for SE events:
        [transcripts_containing_exon1, transcripts_containing_exon2,
        transcripts_containing_exon3]

    Returns
    -------
    isoform1 : set
        Set of transcripts corresponding to isoform1
    isoform2 : set
        Set of transcripts corresponding to isoform2

    Raises
    ------
    """
    isoform1s, isoform2s = None, None
    if splice_type == 'SE':
        isoform1s = set(transcripts[0]) & set(transcripts[2])
        isoform2s = set(transcripts[1]) & isoform1s
    elif splice_type == 'MXE':
        # Isoform 1 is inclusion of the far, second alternative exon
        isoform1s = set(transcripts[0]) & set(transcripts[2]) & \
                    set(transcripts[3])
        # Isoform 2 is inclusion of the near, first alternative exon
        isoform2s = set(transcripts[0]) & set(transcripts[1]) & \
                    set(transcripts[3])
    return isoform1s, isoform2s


def seq_name_to_exon_id(seq_name):
    chr_start_stop, strand = seq_name.split('(')
    chrom, startstop = chr_start_stop.split(':')
    start, stop = startstop.split('-')
    start = int(start)+1
    chr_start_stop = '{}:{}-{}'.format(chrom, start, stop)
    strand = strand.rstrip(')')
    exon = 'exon:{}:{}'.format(chr_start_stop, strand)
    return exon


def isoform_sequences(splice_type, miso_ids):
    """Get the mRNA sequences of the splicing events

    Parameters
    ----------
    splice_type : 'SE' | 'MXE' | 'A5SS' | 'A3SS' | 'RI'
        Name of the type of splicing event
    miso_ids : list-like
        List of miso ID strings, where each exon is delimited by @,
        e.g. for SE: chr1:100:200:+@chr1:300:350:+@chr1:400:500:+


    Returns
    -------


    Raises
    ------
    """
    pass


def isoform_translations(splice_type, miso_ids, gffdb):
    """Get the protein translations (when possible) of the splicing events

     Uses the gffdb to check whether the exons from the miso ID correspond to
     a coding sequence (CDS)

    Parameters
    ----------
    splice_type : 'SE' | 'MXE' | 'A5SS' | 'A3SS' | 'RI'
        Name of the type of splicing event
    miso_ids : list-like
        List of miso ID strings, where each exon is delimited by @,
        e.g. for SE: chr1:100:200:+@chr1:300:350:+@chr1:400:500:+
    gffdb : gffutils.FeatureDB
        A gffutils database, which must have been created by
        create_gffutils_db.create_db because it requires that coding sequences
        aka CDS's are stored like CDS:chr1:100-200:2 where the last 2 is the
        frame of the CDS

    Returns
    -------


    Raises
    ------

    """
    isoform1_seqs = []
    isoform2_seqs = []
    isoform1_translations = defaultdict(list)
    isoform2_translations = defaultdict(list)


    for miso_id in miso_ids:
        exon_ids = miso_id_to_exon_ids(miso_id)
#         print seq1.name

        try:
#             print seq1.name, exon_id1
            exons = [gffdb[exon_id] for exon_id in exon_ids]
        except gffutils.FeatureNotFoundError:
            continue
#         print exon1
#         print

        cds_ids = map(lambda x: ':'.join(['CDS', x.split('exon:')[1]]),
                      exon_ids)
        transcripts = [gffdb.parents(e, featuretype='transcript')
                       for e in exons]
        transcripts_per_isoform = splice_type_isoforms(splice_type,
                                                       transcripts)
        cds_ids_per_isoform = splice_type_exons(splice_type, cds_ids)

        # Remove all overlapping isoforms
        intersect = transcripts_per_isoform[0].intersect(transcripts_per_isoform[1])
        isoform1s = transcripts_per_isoform[0].difference(intersect)
        isoform2s = transcripts_per_isoform[1].difference(intersect)
        transcripts_per_isoform = isoform1s, isoform2s

        for i, (transcript_isoform, cds_isoform) in enumerate(transcripts_per_isoform,
                                                 cds_ids_per_isoform):
            event_isoform = '{}_isoform{}'.format(miso_id, i+1)
            name = '{}_{}'.format(event_isoform, isoform.id)
            # Check if this is the reversed strand
            reverse = exon_ids[0][-1] == '-'


            for t in transcript_isoform:
                cds = list(gffdb.children(t, featuretype='CDS',
                                          reverse=reverse, order_by='start'))
                cds_in_splice_form = [(i, c) for c in cds if
                                      map(lambda x: x.startswith(c.id),
                                          cds_ids)]

                cds_in_correct_order = True
                for (i, c), (j, d) in zip(cds_in_splice_form,
                                          cds_in_splice_form[1:]):
                    if i+1 != j:
                        cds_in_correct_order = False
                if not cds_in_correct_order:
                    break



        for isoform in isoform1s:




            cds_in_splice_form = [c for c in cds if
                                  c.id.startswith(cds_id1)
                                  or c.id.startswith(cds_id3)]
#             print 'cds_in_splice_form', cds_in_splice_form

            if len(cds_in_splice_form) == 2:
                frame1 = cds_in_splice_form[0].frame
                # frame3 = cds_in_splice_form[1].frame

                seq = seq1.seq + seq2.seq

                seq_translated = seq[int(frame1):].translate()
                if seq_translated in isoform1_translations[event_isoform]:
                    continue

#                 print 'name', name
                result_seq = SeqRecord(seq_translated, id=name, description='')
                isoform1_seqs.append(result_seq)
                isoform1_translations[event_isoform].append(seq_translated)
            else:
                isoform1_translations[event_isoform].append('no translation')

        for isoform in isoform2s:
            event_isoform = '{}_isoform2'.format(miso_id)
            name = '{}_{}'.format(event_isoform, isoform.id)
            reverse = exon_id1[-1] == '-'
            cds = gffdb.children(isoform, featuretype='CDS', reverse=reverse, order_by='start')

            cds_in_splice_form = [c for c in cds if c.id.startswith(cds_id1) or c.id.startswith(cds_id2) or c.id.startswith(cds_id3)]
#             print 'cds_in_splice_form', cds_in_splice_form

            if len(cds_in_splice_form) == 3:
                frame1 = cds_in_splice_form[0].frame
                frame2 = cds_in_splice_form[1].frame
                frame3 = cds_in_splice_form[2].frame

                seq = seq1.seq + seq2.seq + seq3.seq
                seq_translated = seq[int(frame1):].translate()
                if seq_translated in isoform2_translations[event_isoform]:
                    continue

#                 print 'name', name
                result_seq = SeqRecord(seq_translated, id=name, description='')
                isoform2_seqs.append(result_seq)
                isoform2_translations[event_isoform].append(seq_translated)
#                 isoform1_names.add(name)
            else:
                isoform2_translations[event_isoform].append('no translation')

    return isoform1_seqs, isoform2_seqs
