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
        full_match = {'$regex': f'^{gene_prefix}$', '$options': 'i'}
        partial_match = {'$regex': f'^{gene_prefix}.', '$options': 'i'}
        internal_match = {'$regex': f'.{gene_prefix}', '$options': 'i'}
        result = self.genes.aggregate([
            { '$match': {'$or': [
                {'ensembl_id': full_match},
                {'symbol': full_match},
                {'alias': full_match},
                {'refseq': full_match},
            ]} },
            { '$set': { 'matchQuality': 10} },
            { '$unionWith': {
                'coll': 'genes',
                'pipeline': [
                    { '$match': {'$or': [
                        {'ensembl_id': partial_match},
                        {'symbol': partial_match},
                        {'alias': partial_match},
                        {'refseq': partial_match},
                    ]} },
                    { '$set': { 'matchQuality': 5} }             
                ]
            } },
            { '$unionWith': {
                'coll': 'genes',
                'pipeline': [
                    { '$match': {'$or': [
                        {'ensembl_id': internal_match},
                        {'symbol': internal_match},
                        {'alias': internal_match},
                        {'refseq': internal_match},
                    ]} },
                    { '$set': { 'matchQuality': 1} }             
                ]
            } },
            { '$sort': { 'matchQuality': -1 } },
            { '$limit': 100 },
            { '$project': {'_id': 0} },
        ])
        return list(result)

    def find_gene(self, gene_name):
        case_insensitive = {'$regex': f'^{gene_name}$', '$options': 'i'}
        search_expr = {'$or': [
            {'ensembl_id': case_insensitive},
            {'symbol': case_insensitive},
            {'refseq': case_insensitive},
        ]}
        result = list(self.genes.find(search_expr, projection=self.no_id, limit=100))
        if not result:
            result = list(self.genes.find({'alias': case_insensitive}, projection=self.no_id, limit=100))
        return result if result else None

    def get_results_by_gene(self, ensembl_id):
        result = self.results.aggregate([
            { '$match': { 'ensembl_id': ensembl_id.upper() } },
            { '$lookup': { 
                'from': 'viruses', 
                'localField': 'virus', 
                'foreignField': 'species',
                'as': 'virus_info',
            }},
            { '$addFields': {
                'virus_genome': {'$first': '$virus_info.genome_composition'},
                'virus_baltimore': { '$first': '$virus_info.baltimore' },
            }},
            { '$project': {
                '_id': 0,
                'virus_info': 0,
            }},
            { '$sort': { 'adj_p': 1 }},
        ])
        return list(result)
        
    def get_multiple_gene_results(self, ensembl_list):
        """
        Assumes that ensembl_list contains valid (ie, uppercase) IDs
        """
        return list(self.results.find(
            {'ensembl_id': {'$in': ensembl_list}},
            projection=self.no_id)
        )

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
                            'usable_samples': {
                                '$filter': {
                                    'input': '$samples',
                                    'cond': {'$eq': [{'$type': '$$this.missing_reason'}, 'missing']}
                                }
                            },
                            'total_samples': '$samples',
                            'usable_controls': {
                                '$filter': {
                                    'input': '$controls',
                                    'cond': {'$eq': [{'$type': '$$this.missing_reason'}, 'missing']}
                                }
                            },
                            'total_controls': '$controls',
                            'missing_sample_reasons': '$missing_sample_reasons',
                            'missing_control_reasons': '$missing_control_reasons',
                            'cell_type': '$cell_type'
                } ] } }},
            { '$set': {
                    'enough_samples': { '$and': [
                        { '$gte': [{'$size': '$usable_samples'}, 2] },
                        { '$gte': [{'$size': '$usable_controls'}, 2] }
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

    def get_analysis_results(self, study, virus, cell_type, platform, geneSet, num_genes=None):
        initial_match = {
            'study': study,
            'virus': virus,
            'bto_id': cell_type,
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
        pipeline = [
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
        ]
        if num_genes is not None:
            pipeline.append({ '$sort': {'adj_p': 1, '_id': 1} })
            pipeline.append({ '$limit': num_genes })
        pipeline.append({ '$project': { '_id': 0 } })
        return list(self.results.aggregate(pipeline))
    
    def num_studies(self):
        return list(self.geo.aggregate([
            { '$match': {'samples.valid_experiment': True} },
            { '$count': 'num_studies' }
        ]))[0]['num_studies']
    
    def num_viruses(self):
        return list(self.geo.aggregate([
            { '$match': {'samples.valid_experiment': True} },
            { '$unwind': {'path': '$samples'} },
            { '$match': {
                'samples.valid_experiment': True,
                'samples.normalized_virus': { '$ne': 'Uninfected' }
            } },
            { '$group': {'_id': '$samples.normalized_virus'} },
            { '$count': 'num_viruses' }
        ]))[0]['num_viruses']