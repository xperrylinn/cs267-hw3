#pragma once
#include <map>
#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>

struct HashMap {

    using dobj_kmer_map_t = upcxx::dist_object<std::vector<kmer_pair>>;
    dobj_kmer_map_t data;

    using dobj_map_usage_t = upcxx::dist_object<std::vector<int>>;
    dobj_map_usage_t used;

    upcxx::global_ptr<int> my_size_gptr;

    // Original code
    // std::vector<kmer_pair> data;
    // std::vector<int> used;

    // size_t my_size;

    size_t size() const noexcept;

    HashMap(size_t size);

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    upcxx::future<bool> insert(const kmer_pair& kmer);
    upcxx::future<bool> find(const pkmer_t& key_kmer, kmer_pair& val_kmer);

    // Helper functions
    uint64_t get_target_rank(const kmer_pair& kmer);

    // Write and read to a logical data slot in the table.
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);
};

HashMap::HashMap(size_t size) {
    my_size_gptr = upcxx::new_<int>(size);
    data->resize(size);
    used->resize(size, 0);
}

upcxx::future<bool> HashMap::insert(const kmer_pair& kmer) {
    upcxx::future<bool> fut_result = upcxx::rpc(
        get_target_rank(kmer),
        [](dobj_kmer_map_t &data, dobj_map_usage_t &used, upcxx::global_ptr<int>& my_size, const kmer_pair& kmer) -> bool {
            uint64_t hash = kmer.hash();
            uint64_t probe = 0;
            bool success = false;
            do {
                uint64_t slot = (hash + probe++) % *my_size.local();
                // success = request_slot(slot);
                if (used->data()[slot] != 0) {
                    success = false;
                } else {
                    used->data()[slot] = 1;
                    success = true;
                }
                if (success) {
                    // write_slot(slot, kmer);
                    data->data()[slot] = kmer;
                }
            } while (!success && probe < *my_size.local());
            return success;
        },
        data,
        used,
        my_size_gptr,
        kmer
    );
    return fut_result;
}

upcxx::future<bool> HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    upcxx::future<bool> fut_result = upcxx::rpc(
        get_target_rank(val_kmer),
        [](dobj_kmer_map_t &data, dobj_map_usage_t &used, upcxx::global_ptr<int>& my_size, const pkmer_t& key_kmer, kmer_pair& val_kmer) -> bool {
            uint64_t hash = key_kmer.hash();
            uint64_t probe = 0;
            bool success = false;
            do {
                uint64_t slot = (hash + probe++) % *my_size.local();
                bool is_slot_used = used->data()[slot] != 0;
                // if (slot_used(slot)) {
                if (is_slot_used) {
                    // val_kmer = read_slot(slot);
                    val_kmer = data->data()[slot];
                    if (val_kmer.kmer == key_kmer) {
                        success = true;
                    }
                }
            } while (!success && probe < *my_size.local());
            return success;
        },
        data,
        used,
        my_size_gptr,
        key_kmer,
        val_kmer
    );
    return fut_result;
}

uint64_t get_target(const kmer_pair& kmer) {
    return kmer.hash() % upcxx::rank_n();
}

bool HashMap::slot_used(uint64_t slot) { return used->data()[slot] != 0; }

void HashMap::write_slot(uint64_t slot, const kmer_pair& kmer) { data->data()[slot] = kmer; }

kmer_pair HashMap::read_slot(uint64_t slot) { return data->data()[slot]; }

bool HashMap::request_slot(uint64_t slot) {
    if (used->data()[slot] != 0) {
        return false;
    } else {
        used->data()[slot] = 1;
        return true;
    }
}

size_t HashMap::size() const noexcept { return my_size; }
