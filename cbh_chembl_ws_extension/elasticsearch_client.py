
from django.conf import settings
import elasticsearch

import time
try:
    ES_PREFIX = settings.ES_PREFIX
except AttributeError:
    ES_PREFIX = "dev"

def get_temp_index_name(request, multi_batch_id):
    return "%s__temp_multi_batch__%s__%s" % (ES_PREFIX, request.session.session_key, str(multi_batch_id))



def get(index_name, es_request_body, bundledata):
    es = elasticsearch.Elasticsearch()
    result = es.search(index_name, body=es_request_body)
    data = [res["_source"] for res in result["hits"]["hits"]]
    
    bundledata["meta"] = {"totalCount" : result["hits"]["total"]}
    bundledata["objects"] = data
    return bundledata



def create_temporary_index(batches, multi_batch_id, request):
    es = elasticsearch.Elasticsearch()
    t = time.time()
    store_type = "memory"
    if len(batches) > 100000:
        store_type = "niofs"
    create_body = {
        "settings": {
            "index.store.type": store_type
        },
         "mappings" : {
            "_default_" : {
               "_all" : {"enabled" : False},
               "dynamic_templates" : [ {
                 "string_fields" : {
                   "match" : "ctab|std_ctab|canonical_smiles|original_smiles",
                   "match_mapping_type" : "string",
                   "mapping" : {
                        "type" : "string", "store" : "no", "include_in_all" : False        
                   }
                 }
               } ,


                {
                 "string_fields" : {
                   "match" : "*",
                   "match_mapping_type" : "string",
                   "mapping" : {
                     "type" : "string","store" : "no", "index_options": "docs","index" : "analyzed", "omit_norms" : True,
                       "fields" : {
                         "raw" : {"type": "string","store" : "no", "index" : "not_analyzed", "ignore_above" : 256}
                       }
                   }
                 }
               } 
            ]
        }
        }
    }
    index_name = get_temp_index_name(request, multi_batch_id)
    
    es.indices.create(
            index_name,
            body=create_body,
            ignore=400)
    
    bulk_items = []
    for item in batches:
        bulk_items.append({
                            "index" :
                                {
                                    "_id": str(item["id"]), 
                                    "_index": index_name,
                                    "_type": "%d" % multi_batch_id
                                }
                            })
        bulk_items.append(item)
    #Data is not refreshed!
    es.bulk(body=bulk_items)
