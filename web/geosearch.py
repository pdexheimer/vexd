from pymongo import MongoClient


class GeoSearch:

    def __init__(self, host='localhost', port=27017, user=None, password=None):
        self.client = MongoClient(
            host=host,
            port=port,
            username=user,
            password=password,
            connect=False
        )
        self.geo = self.client.vexd.geo
        self.virus = self.client.vexd.viruses
        self.genes = self.client.vexd.genes
        self.results = self.client.vexd.results
        self.bto = self.client.vexd.bto
        self.no_id = {'_id': False}

    def get_gse_info(self, gse_id):
        return self.geo.find_one({'id': gse_id.upper()}, projection=self.no_id)

    def get_gene_info(self, ensembl_id):
        return self.genes.find_one({'ensembl_id': ensembl_id.upper()}, projection=self.no_id)

    def search_gse(self, virus_name=None, platform=None, with_pubmed=False,
                   only_valid=True, human_only_virus=True):
        search_expr = {}
        # If multiple clauses are requested in the samples array, use $elemMatch
        # to ensure that (at least) one sample meets all the conditions
        if virus_name is not None and virus_name:
            if not self.is_virus_name_canonical(virus_name):
                virus_name = self.resolve_virus_alias(
                    virus_name, human_only=human_only_virus)
            else:
                virus_name = [virus_name, ]
            if only_valid:
                search_expr['samples'] = {
                    '$elemMatch': {
                        'normalized_virus': {'$in': virus_name},
                        'valid_experiment': True
                    }
                }
            else:
                search_expr['samples.normalized_virus'] = {'$in': virus_name}
        elif only_valid:
            search_expr['samples.valid_experiment'] = True
        if platform is not None:
            search_expr['platform.type'] = platform
        if with_pubmed:
            search_expr['pubmed'] = {'$exists': True}
        return list(self.geo.find(search_expr, projection=self.no_id))

    def count_by_virus(self):
        result = self.geo.aggregate([
            {'$unwind': '$samples'},
            {'$match': {'samples.valid_experiment': True}},
            {'$group': {
                '_id': '$samples.normalized_virus',
                'sample_count': {'$sum': 1},
                'studies': {'$addToSet': '$id'}
            }},
            {'$set': {'study_count': {'$size': '$studies'}}},
            {'$sort': {
                'study_count': -1,
                'sample_count': -1
            }},
        ])
        return list(result)

    def resolve_virus_alias(self, virus_alias, human_only=True):
        case_insensitive = {'$regex': f'^{virus_alias}$', '$options': 'i'}
        search_expr = {'$or': [
            {'species': case_insensitive},
            {'name': case_insensitive},
            {'aliases': case_insensitive},
        ]}
        if human_only:
            search_expr['human_infecting'] = True
        return [x['species'] for x in self.virus.find(search_expr, projection={'species': True})]

    def is_virus_name_canonical(self, virus_name):
        if virus_name == 'Uninfected':
            return True
        return self.get_virus(virus_name) is not None

    def get_virus(self, virus_name):
        return self.virus.find_one({'species': virus_name}, projection=self.no_id)

    def find_virus_by_prefix(self, virus_prefix, human_only=True):
        case_insensitive = {'$regex': f'^{virus_prefix}', '$options': 'i'}
        search_expr = {'$or': [
            {'species': case_insensitive},
            {'name': case_insensitive},
            {'aliases': case_insensitive},
        ]}
        if human_only:
            search_expr['human_infecting'] = True
        return list(self.virus.find(search_expr, projection=self.no_id))

    def find_gene_by_prefix(self, gene_prefix):
        case_insensitive = {'$regex': f'^{gene_prefix}', '$options': 'i'}
        search_expr = {'$or': [
            {'ensembl_id': case_insensitive},
            {'symbol': case_insensitive},
            {'alias': case_insensitive},
            {'refseq': case_insensitive},
        ]}
        return list(self.genes.find(search_expr, projection=self.no_id, limit=100))

    def get_results_by_gene(self, ensembl_id):
        return list(self.results.find({'ensembl_id': ensembl_id.upper()}, sort=[('adj_p', 1)]))

    def get_valid_comparisons_by_virus(self, virus_name):
        result = self.results.aggregate([
            {'$match': {'virus': virus_name}},
            {'$group': {
                '_id': {
                    'study': '$study',
                    'cell_type': '$cell_type',
                    'platform': '$platform'
                }
            }},
            {'$replaceRoot': {'newRoot': '$_id'}}
        ])
        return list(result)

    def get_gene_results_by_virus(self, virus_name, pval_cutoff=0.05, logfc_cutoff=1, num_genes=10):
        # This is a truly massive aggregation.  The stages are:
        # 1) match - select only result for the virus
        # 2) group - figure out the unique combos of study/cell_type/platform
        # 3) lookup - get the study title
        # 4) replaceRoot - promote the results from the stage 3 lookup to the top level
        # 5) lookup - count the number of significant genes
        # 6) set - promote the (single) object of counts out of the array
        # 7) lookup - get the top num_genes upregulated genes
        # 8) lookup - get the top num_genes downregulated genes
        # 9) sort - by decreasing number of significant genes
        result = self.results.aggregate([
            {'$match': {'virus': virus_name}},
            {'$group': {
                '_id': {
                    'study': '$study',
                    'cell_type': '$cell_type',
                    'platform': '$platform',
                    'virus': '$virus',
                    'bto_id': '$bto_id',
                    'bto_name': '$bto_name',
                    'cellosaurus_id': '$cellosaurus_id',
                    'cellosaurus_name': '$cellosaurus_name'
                },
                'infect_count': {'$first': '$infect_count'},
                'control_count': {'$first': '$control_count'},
            }
            },
            {'$lookup': {
                'from': 'geo',
                'let': {'studyid': '$_id.study'},
                'pipeline': [
                    {'$match': {'$expr': {'$eq': ['$id', '$$studyid']}}},
                    {'$project': {'_id': 0, 'title': 1}}
                ],
                'as': 'geo'
            }
            }, {
                '$replaceRoot': {
                    'newRoot': {'$mergeObjects': [
                        {'$arrayElemAt': ['$geo', 0]},
                        '$_id',
                        {'infect_count': '$infect_count',
                            'control_count': '$control_count'}
                    ]}
                }
            }, {
                '$lookup': {
                    'from': 'results',
                    'let': {
                        'study': '$study',
                        'cell_type': '$cell_type',
                        'platform': '$platform'
                    },
                    'pipeline': [
                        {
                            '$match': {
                                '$expr': {
                                    '$and': [
                                        {'$eq': ['$$study', '$study']},
                                        {'$eq': ['$$cell_type', '$cell_type']},
                                        {'$eq': ['$$platform', '$platform']},
                                        {'$lte': ['$adj_p', pval_cutoff]},
                                        {'$or': [
                                            {'$gte': ['$logfc', logfc_cutoff]},
                                            {'$lte': ['$logfc', -logfc_cutoff]}
                                        ]}
                                    ]
                                }
                            }
                        }, {
                            '$group': {
                                '_id': None,
                                'signif': {'$sum': 1},
                                'down': {'$sum': {'$cond': [{'$lt': ['$logfc', 0]}, 1, 0]}},
                                'up': {'$sum': {'$cond': [{'$gt': ['$logfc', 0]}, 1, 0]}}
                            }
                        },
                        {'$project': {'_id': 0}}
                    ],
                    'as': 'counts'
                }
            },
            {'$set': {'counts': {'$mergeObjects': [
                {'signif': 0, 'up': 0, 'down': 0},
                {'$arrayElemAt': ['$counts', 0]}
            ]}}},
            {
                '$lookup': {
                    'from': 'results',
                    'let': {
                        'study': '$study',
                        'cell_type': '$cell_type',
                        'platform': '$platform'
                    },
                    'pipeline': [
                        {
                            '$match': {
                                '$expr': {
                                    '$and': [
                                        {'$eq': ['$$study', '$study']},
                                        {'$eq': ['$$cell_type', '$cell_type']},
                                        {'$eq': ['$$platform', '$platform']},
                                        {'$lte': ['$adj_p', pval_cutoff]},
                                        {'$gte': ['$logfc', logfc_cutoff]}
                                    ]
                                }
                            }
                        },
                        {'$sort': {'logfc': -1, 'adj_p': 1}},
                        {'$limit': num_genes},
                        {
                            '$lookup': {
                                'from': 'genes',
                                'localField': 'ensembl_id',
                                'foreignField': 'ensembl_id',
                                'as': 'symbol'
                            }
                        }, {
                            '$project': {
                                'ensembl_id': 1,
                                'adj_p': 1,
                                'logfc': 1,
                                'fc': 1,
                                '_id': 0,
                                'symbol': {'$arrayElemAt': ['$symbol', 0]}
                            }
                        }, {
                            '$addFields': {
                                'description': '$symbol.description',
                                'symbol': '$symbol.symbol'
                            }
                        }
                    ],
                    'as': 'up'
                }
            },
            {
                '$lookup': {
                    'from': 'results',
                    'let': {
                        'study': '$study',
                        'cell_type': '$cell_type',
                        'platform': '$platform'
                    },
                    'pipeline': [
                        {
                            '$match': {
                                '$expr': {
                                    '$and': [
                                        {'$eq': ['$$study', '$study']},
                                        {'$eq': ['$$cell_type', '$cell_type']},
                                        {'$eq': ['$$platform', '$platform']},
                                        {'$lte': ['$adj_p', pval_cutoff]},
                                        {'$lte': ['$logfc', -logfc_cutoff]}
                                    ]
                                }
                            }
                        },
                        {'$sort': {'logfc': 1, 'adj_p': 1}},
                        {'$limit': num_genes},
                        {
                            '$lookup': {
                                'from': 'genes',
                                'localField': 'ensembl_id',
                                'foreignField': 'ensembl_id',
                                'as': 'symbol'
                            }
                        }, {
                            '$project': {
                                'ensembl_id': 1,
                                'adj_p': 1,
                                'logfc': 1,
                                'fc': 1,
                                '_id': 0,
                                'symbol': {'$arrayElemAt': ['$symbol', 0]}
                            }
                        }, {
                            '$addFields': {
                                'description': '$symbol.description',
                                'symbol': '$symbol.symbol'
                            }
                        }
                    ],
                    'as': 'down'
                }
            },
            {'$sort': {'counts.signif': -1}}
        ])
        return list(result)

    def get_significant_genes_by_virus(self, virus_name, pval_cutoff=0.05, logfc_cutoff=1):
        # Another aggregation:
        # 1) match - select the correct virus
        # 2) group - by probeset id.  Retain the ensembl id as a separate element
        # 3) match - select only the probesets that are significant in at least one study
        # 4) lookup - the gene symbol
        # 5) set - extract the looked up gene info from the array (and set a default)
        # 6) set - extract the gene symbol from the gene info
        # 7) unwind - now that we only have sig genes, unwind them for sorting
        # 8) sort - to get the studies array into a consistent order
        # 9) group - everything by probeset again
        # 10) sort - the results by gene symbol
        result = self.results.aggregate([
            {'$match': {'virus': virus_name}},
            {'$group': {
                '_id': '$probeset_id',
                'ensembl_id': {'$first': '$ensembl_id'},
                'studies': {'$push': '$$ROOT'}
            }},
            {'$match': {
                'studies': {
                    '$elemMatch': {
                        'adj_p': {'$lte': pval_cutoff},
                        '$or': [
                            {'logfc': {'$gte': logfc_cutoff}},
                            {'logfc': {'$lte': -logfc_cutoff}}
                        ]
                    }
                }
            }},
            {'$lookup': {
                'from': 'genes',
                'localField': 'ensembl_id',
                'foreignField': 'ensembl_id',
                'as': 'symbol'
            }},
            {'$set': {
                'symbol': {
                    '$mergeObjects': [
                        {'symbol': ''},
                        {'$arrayElemAt': ['$symbol', 0]}
                    ]
                }
            }},
            {'$set': {'symbol': '$symbol.symbol'}},
            {'$unwind': {'path': '$studies'}},
            {'$sort': {
                'studies.study': 1,
                'studies.cell_type': 1,
                'studies.platform': 1
            }},
            {'$group': {
                '_id': '$_id',
                'ensembl_id': {'$first': '$ensembl_id'},
                'symbol': {'$first': '$symbol'},
                'studies': {'$push': '$studies'}
            }},
            {'$sort': {'symbol': 1}}
        ], allowDiskUse=True)
        return list(result)

    def get_studies_by_virus(self, virus_name):
        # Small aggregation this time:
        # 1) match - select the virus
        # 2) group - get the unique study/cell type/platform tuples
        # 3) project - reformat the results
        # 4) sort - put them in a consistent order
        result = self.results.aggregate([
            {'$match': {'virus': virus_name}},
            {'$group': {'_id': {
                'study': '$study',
                'cell_type': '$cell_type',
                'bto_name': '$bto_name',
                'platform': '$platform'
            }}},
            {'$project': {
                'study': '$_id.study',
                'cell_type': '$_id.cell_type',
                'bto_name': '$_id.bto_name',
                'platform': '$_id.platform',
                '_id': 0
            }},
            {'$sort': {
                'study': 1,
                'cell_type': 1,
                'platform': 1
            }}
        ])
        return list(result)

    def search_studies(self, virus_name, tissue, use_descendants, pval_cutoff=0.05, logfc_cutoff=1):
        #  1) match: Restrict to only studies of interest
        #            Partially redundant with #3, but will make #2 faster
        #  2) unwind: Work on the sample level
        #  3) match: Restrict to only samples of interest
        #  4) group: By study, cell type, platform, and virus
        #  5) lookup: Collect the controls for this comparison
        #  6) lookup: Count the number of significant genes
        #  7) lookup: Get the study title
        #  8) lookup: Collect the missing_sample_reasons
        #  9) lookup: Collect the missing_control_reasons
        # 10) replaceRoot: Promote fields from #4-9 (except 5) to top level,
        #                  remove samples/controls that have a 'missing_reason' field
        # 11) set: Flag the comparisons which don't have enough samples
        # 12) sort: by # significant genes (descending), move the small comparisons to the end
        initial_filter = {'samples.valid_experiment': True}
        if virus_name:
            initial_filter['samples.normalized_virus'] = virus_name
        if tissue:
            if use_descendants:
                target_ids = self.get_descendants_of_bto_term(tissue)
                initial_filter['samples.bto_id'] = { '$in': target_ids }
            else:
                initial_filter['samples.bto_id'] = tissue
        second_filter = dict(initial_filter)
        if not virus_name:
            second_filter['samples.normalized_virus'] = { '$ne': 'Uninfected' }
        result = self.geo.aggregate([
            { '$match': initial_filter },
            { '$unwind': {'path': '$samples'} },
            { '$match': second_filter },
            { '$group': {
                    '_id': {
                        'study': '$id',
                        'bto_name': '$samples.bto_name',
                        'virus': '$samples.normalized_virus',
                        'platform': '$samples.cdf',
                        'bto_id': '$samples.bto_id',
                    },
                    'samples': { '$push': '$samples' },
                    'cell_type': { '$first': '$samples.cell_type' }
                } },
            { '$lookup': {
                    'from': 'geo',
                    'as': 'controls',
                    'let': {
                        'study': '$_id.study',
                        'cell_type': '$_id.bto_name',
                        'platform': '$_id.platform'
                    },
                    'pipeline': [
                        { '$match': { '$expr': {'$eq': ['$id', '$$study'] } } },
                        { '$unwind': {'path': '$samples'} }, 
                        { '$match': { '$expr': {
                            '$and': [
                                { '$eq': ['$samples.normalized_virus', 'Uninfected'] },
                                { '$eq': ['$samples.valid_experiment', True] },
                                { '$eq': ['$id', '$$study'] },
                                { '$eq': ['$samples.cdf', '$$platform'] },
                                { '$eq': ['$samples.bto_name', '$$cell_type'] }
                            ] } } }, 
                        { '$replaceRoot': {'newRoot': '$samples'} }
                    ]
                } }, 
            { '$lookup': {
                    'from': 'results',
                    'let': {
                        'study': '$_id.study',
                        'cell_type': '$_id.bto_name',
                        'platform': '$_id.platform',
                        'virus': '$_id.virus'
                    },
                    'pipeline': [
                        { '$match': { '$expr': {
                            '$and': [
                                { '$eq': ['$$study', '$study'] },
                                { '$eq': ['$$cell_type', '$bto_name'] },
                                { '$eq': ['$$platform', '$platform'] },
                                { '$eq': ['$$virus', '$virus'] },
                                { '$lte': ['$adj_p', pval_cutoff] },
                                { '$or': [
                                    { '$gte': ['$logfc', logfc_cutoff] },
                                    { '$lte': ['$logfc', -logfc_cutoff] }
                                ] }
                            ] } } },
                        { '$group': {
                            '_id': None,
                            'signif': { '$sum': 1 },
                            'down': { '$sum': {'$cond': [{'$lt': ['$logfc', 0]}, 1, 0]} },
                            'up': { '$sum': {'$cond': [{'$gt': ['$logfc', 0]}, 1, 0]} }
                            } },
                        { '$project': {'_id': 0} }
                    ],
                    'as': 'gene_count'
                } },
            { '$lookup': {
                    'from': 'geo',
                    'let': {'studyid': '$_id.study'},
                    'pipeline': [
                        {'$match': {'$expr': {'$eq': ['$id', '$$studyid']}}},
                        {'$project': {'_id': 0, 'title': 1}}
                    ],
                    'as': 'geo' 
                } },
            { '$lookup': {
                    'from': 'geo',
                    'let': {
                        'study': '$_id.study',
                        'cell_type': '$_id.bto_name',
                        'platform': '$_id.platform',
                        'virus': '$_id.virus'
                        },
                    'pipeline': [
                        {'$match': { '$expr': {'$eq': ['$id', '$$study'] }}},
                        {'$unwind': { 'path': '$samples' }},
                        {'$match': {'$expr': {
                            '$and': [
                                {'$eq': ['$samples.valid_experiment', True]},
                                {'$eq': ['$samples.normalized_virus', '$$virus']},
                                {'$eq': ['$samples.bto_name', '$$cell_type']},
                                {'$eq': ['$samples.cdf', '$$platform']}
                            ]
                        } } },
                        {'$group': {
                            '_id': '$samples.missing_reason',
                            'count': {'$sum': 1}
                        }}
                    ],
                    'as': 'missing_sample_reasons'
                } },
            { '$lookup': {
                    'from': 'geo',
                    'let': {
                        'study': '$_id.study',
                        'cell_type': '$_id.bto_name',
                        'platform': '$_id.platform'
                        },
                    'pipeline': [
                        {'$match': { '$expr': {'$eq': ['$id', '$$study'] }}},
                        {'$unwind': { 'path': '$samples' }},
                        {'$match': {'$expr': {
                            '$and': [
                                {'$eq': ['$samples.valid_experiment', True]},
                                {'$eq': ['$samples.normalized_virus', 'Uninfected']},
                                {'$eq': ['$samples.bto_name', '$$cell_type']},
                                {'$eq': ['$samples.cdf', '$$platform']}
                            ]
                        } } },
                        {'$group': {
                            '_id': '$samples.missing_reason',
                            'count': {'$sum': 1}
                        }}
                    ],
                    'as': 'missing_control_reasons'
                } },
            { '$replaceRoot': {'newRoot': {
                    '$mergeObjects': [
                        '$_id', 
                        {'$arrayElemAt': ['$gene_count', 0]},
                        {'$arrayElemAt': ['$geo', 0]},
                        {
                            'pf_samples': {
                                '$filter': {
                                    'input': '$samples',
                                    'cond': {'$eq': [{'$type': '$$this.missing_reason'}, 'missing']}
                                }
                            },
                            'samples': '$samples',
                            'pf_controls': {
                                '$filter': {
                                    'input': '$controls',
                                    'cond': {'$eq': [{'$type': '$$this.missing_reason'}, 'missing']}
                                }
                            },
                            'controls': '$controls',
                            'missing_sample_reasons': '$missing_sample_reasons',
                            'missing_control_reasons': '$missing_control_reasons',
                            'cell_type': '$cell_type'
                } ] } }},
            { '$set': {
                    'enough_samples': { '$and': [
                        { '$gte': [{'$size': '$pf_samples'}, 2] },
                        { '$gte': [{'$size': '$pf_controls'}, 2] }
                    ] },
                    'study_order': { '$toInt': { '$ltrim': {
                        'input': '$study',
                        'chars': 'GSE'
                    } } }
                }},
            { '$sort': { 'signif': -1, 'enough_samples': -1, 'study_order': 1 } }
        ])
        return list(result)
    
    def get_descendants_of_bto_term(self, bto_id):
        result = self.bto.aggregate([
            { '$match': {'id': bto_id} },
            { '$graphLookup': {
                'from': 'bto',
                'startWith': '$id',
                'connectFromField': 'descendants',
                'connectToField': 'id',
                'as': 'kids'
            }},
            { '$project': {
                '_id': 0,
                'results': '$kids.id'
            }}
        ])
        return result.next()['results']

    def get_analysis_results(self, study, virus, cell_type, platform, geneSet):
        initial_match = {
            'study': study,
            'virus': virus,
            'bto_name': cell_type,
        }
        if platform:
            initial_match['platform'] = platform
        if geneSet != 'all':
            initial_match['adj_p'] = {'$lte': 0.05}
        if geneSet == 'sig':
            initial_match['$or'] = [
                {'logfc': {'$gte': 1}},
                {'logfc': {'$lte': -1}}
            ]
        elif geneSet == 'up':
            initial_match['logfc'] = {'$gte': 1}
        elif geneSet == 'down':
            initial_match['logfc'] = {'$lte': -1}
        result = self.results.aggregate([
            { '$match': initial_match },
            { '$lookup': {
                'from': 'genes',
                'localField': 'ensembl_id',
                'foreignField': 'ensembl_id',
                'as': 'symbol'
            }},
            { '$set': { 'symbol': {'$mergeObjects': [
                { 'symbol': '' },
                { '$arrayElemAt': ['$symbol', 0] }
            ]}}},
            { '$set': { 'symbol': '$symbol.symbol' } },
            { '$project': { '_id': 0 } }
        ])
        return list(result)