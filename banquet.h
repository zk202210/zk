#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "banquet_instances.h"
#include "types.h"

// crypto api
banquet_keypair_t banquet_keygen(const banquet_instance_t &instance);

banquet_signature_t banquet_sign(const banquet_instance_t &instance,
                                 const banquet_keypair_t &keypair,
                                 const uint8_t *message, size_t message_len);

banquet_group_signature_t banquet_group_sign(const banquet_instance_t &instance,
                                 const banquet_keypair_t &keypair,
                                 const std::vector<banquet_keypair_t> &gkey,
                                 const uint8_t *message, size_t message_len);

bool banquet_verify(const banquet_instance_t &instance,
                    const std::vector<uint8_t> &pk,
                    const banquet_signature_t &signature,
                    const uint8_t *message, size_t message_len);

bool banquet_group_verify(const banquet_instance_t &instance,
                          const banquet_group_signature_t &signature,
                          const std::vector<banquet_keypair_t> &gkey,
                          const uint8_t *message, size_t message_len);

std::vector<uint8_t>
banquet_serialize_signature(const banquet_instance_t &instance,
                            const banquet_signature_t &signature);
banquet_signature_t
banquet_deserialize_signature(const banquet_instance_t &instance,
                              const std::vector<uint8_t> &signature);
std::vector<uint8_t>
banquet_serialize_group_signature(const banquet_instance_t &instance,
                                  size_t gsize,
                                  const banquet_group_signature_t &signature);
banquet_group_signature_t
banquet_deserialize_group_signature(const banquet_instance_t &instance,
                                    const std::vector<uint8_t> &signature);
