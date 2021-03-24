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
        self.geo = self.client.geo_test_db.geo
        self.virus = self.client.virus_db.viruses
        self.genes = self.client.geo_test_db.genes
        self.results = self.client.geo_test_db.results
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
                virus_name = self.resolve_virus_alias(virus_name, human_only=human_only_virus)
            else:
                virus_name = [virus_name,]
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
        return [x for x in self.geo.find(search_expr, projection=self.no_id)]
    
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
        return [x for x in result]

    def resolve_virus_alias(self, virus_alias, human_only=True):
        case_insensitive = {'$regex': f'^{virus_alias}$', '$options': 'i'}
        search_expr = {'$or': [
            {'species': case_insensitive},
            {'name': case_insensitive},
            {'aliases': case_insensitive},
        ]}
        if human_only:
            search_expr['human_infecting'] = True
        return [ x['species'] for x in self.virus.find(search_expr, projection={'species': True}) ]

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
        return [ x for x in self.virus.find(search_expr, projection=self.no_id) ]

    def find_gene_by_prefix(self, gene_prefix):
        case_insensitive = {'$regex': f'^{gene_prefix}', '$options': 'i'}
        search_expr = {'$or': [
            {'ensembl_id': case_insensitive},
            {'symbol': case_insensitive},
            {'alias': case_insensitive},
            {'refseq': case_insensitive},
        ]}
        return [ x for x in self.genes.find(search_expr, projection=self.no_id, limit=100) ]
    
    def get_results_by_gene(self, ensembl_id):
        return [ x for x in self.results.find({'ensembl_id': ensembl_id}, sort=[('adj_p', 1)]) ]