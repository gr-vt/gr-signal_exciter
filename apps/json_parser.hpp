/* -*- c++ -*- */
/*
 * Copyright 2016 Bill Clark.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_SIGNAL_EXCITER_JSON_PARSER_HPP
#define INCLUDED_SIGNAL_EXCITER_JSON_PARSER_HPP

#include <cstring>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include <jsoncpp/json/json.h>
#include <signal_exciter/random_signal_config.h>

namespace signal_exciter
{
namespace json
{

  class Signal_Parameters
  {
   private:
    bool                                d_Loaded;
    gr::signal_exciter::sig_params      d_Sig_Params;

    std::vector<char> d_load_check;
    void load_check(){
      bool loaded = true;
      size_t idx = 0;
      while(loaded && (idx < 3)){
        loaded = loaded && d_load_check[idx++];
      }
      if(loaded && (idx>=3)){ d_Loaded = true; }
    }

    void parse_cwmorse(Json::Value& signal){
      d_Sig_Params.char_per_word = signal.get("Characters per Word", 4).asInt();
      d_Sig_Params.words_per_minute = signal.get("Words per Minute", 120).asFloat();
      d_Sig_Params.base_word = signal.get("Base Word", false).asBool();
      //Base Word: (PARIS : false; CODEX : true)
      d_Sig_Params.sps = signal.get("Interpolation", 1).asInt();
      d_Sig_Params.pulse_len = signal.get("Interpolation Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Interpolation Taps"].begin();
          (it != signal["Interpolation Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      float gen_at =
        (d_Sig_Params.words_per_minute/(d_Sig_Params.base_word ? 1. : 1.2)) *
        d_Sig_Params.sps;
      float sanity = (d_Sig_Params.fs - gen_at)/gen_at;
      if(!((sanity > -1.19e-7)&&(sanity < 1.19e-7))){
        std::cerr << "signal_exciter:: json_parser:: signal_parameters:: " <<
            "parse_cwmorse:: (WARN) 'Generation Sample Rate' != 'parameter space' " <<
            "defaulting to 'Generation Sample Rate' <= 'parameter space'" << std::endl;
        d_Sig_Params.fs = gen_at;
      }
      d_load_check[1] = 1;
    }

    void parse_dsb(Json::Value& signal){
      d_Sig_Params.am_norm = true;
      d_Sig_Params.mod_idx = signal.get("Modulation Index", 0.5).asFloat();
      d_Sig_Params.components = signal.get("Spectral Shaping Components", 0).asUInt();
      d_Sig_Params.spectral_len = signal.get("Spectral Shaping Tap Length", 256).asUInt();
      d_Sig_Params.f_max = signal.get("Spectral Shaping Max Frequency", 48e3).asFloat();
      d_Sig_Params.mu = std::vector<float>(d_Sig_Params.components,0.);
      d_Sig_Params.sigma = std::vector<float>(d_Sig_Params.components,0.);
      d_Sig_Params.weight = std::vector<float>(d_Sig_Params.components,0.);
      size_t idx = 0;
      for(Json::Value::iterator it1 = signal["Spectral Shaping Means"].begin(),
          it2 = signal["Spectral Shaping Std"].begin(),
          it3 = signal["Spectral Shaping Weights"].begin();
          ((it1 != signal["Spectral Shaping Means"].end()) &&
          (it2 != signal["Spectral Shaping Std"].end()) &&
          (it3 != signal["Spectral Shaping Weights"].end())) &&
          (idx < d_Sig_Params.components);
          it1++, it2++, it3++, idx++){
        d_Sig_Params.mu[idx] = (*it1).asFloat();
        d_Sig_Params.sigma[idx] = (*it2).asFloat();
        d_Sig_Params.weight[idx] = (*it3).asFloat();
      }
      d_Sig_Params.sps = signal.get("Interpolation", 1).asInt();
      d_Sig_Params.pulse_len = signal.get("Interpolation Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      idx = 0;
      for(Json::Value::iterator it = signal["Interpolation Taps"].begin();
          (it != signal["Interpolation Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      float gen_at = d_Sig_Params.f_max*d_Sig_Params.sps;
      float sanity = (d_Sig_Params.fs - gen_at)/gen_at;
      if(!((sanity > -1.19e-7)&&(sanity < 1.19e-7))){
        std::cerr << "signal_exciter:: json_parser:: signal_parameters:: " <<
            "parse_dsb:: (WARN) 'Generation Sample Rate' != 'parameter space' " <<
            "defaulting to 'Generation Sample Rate' <= 'parameter space'" << std::endl;
        d_Sig_Params.fs = gen_at;
      }
      d_load_check[1] = 1;
    }

    void parse_dsbsc(Json::Value& signal){
      d_Sig_Params.am_norm = true;
      d_Sig_Params.mod_idx = signal.get("Modulation Index", 0.5).asFloat();
      d_Sig_Params.components = signal.get("Spectral Shaping Components", 0).asUInt();
      d_Sig_Params.spectral_len = signal.get("Spectral Shaping Tap Length", 256).asUInt();
      d_Sig_Params.f_max = signal.get("Spectral Shaping Max Frequency", 48e3).asFloat();
      d_Sig_Params.mu = std::vector<float>(d_Sig_Params.components,0.);
      d_Sig_Params.sigma = std::vector<float>(d_Sig_Params.components,0.);
      d_Sig_Params.weight = std::vector<float>(d_Sig_Params.components,0.);
      size_t idx = 0;
      for(Json::Value::iterator it1 = signal["Spectral Shaping Means"].begin(),
          it2 = signal["Spectral Shaping Std"].begin(),
          it3 = signal["Spectral Shaping Weights"].begin();
          ((it1 != signal["Spectral Shaping Means"].end()) &&
          (it2 != signal["Spectral Shaping Std"].end()) &&
          (it3 != signal["Spectral Shaping Weights"].end())) &&
          (idx < d_Sig_Params.components);
          it1++, it2++, it3++, idx++){
        d_Sig_Params.mu[idx] = (*it1).asFloat();
        d_Sig_Params.sigma[idx] = (*it2).asFloat();
        d_Sig_Params.weight[idx] = (*it3).asFloat();
      }
      d_Sig_Params.sps = signal.get("Interpolation", 1).asInt();
      d_Sig_Params.pulse_len = signal.get("Interpolation Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      idx = 0;
      for(Json::Value::iterator it = signal["Interpolation Taps"].begin();
          (it != signal["Interpolation Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      float gen_at = d_Sig_Params.f_max*d_Sig_Params.sps;
      float sanity = (d_Sig_Params.fs - gen_at)/gen_at;
      if(!((sanity > -1.19e-7)&&(sanity < 1.19e-7))){
        std::cerr << "signal_exciter:: json_parser:: signal_parameters:: " <<
            "parse_dsbsc:: (WARN) 'Generation Sample Rate' != 'parameter space' " <<
            "defaulting to 'Generation Sample Rate' <= 'parameter space'" << std::endl;
        d_Sig_Params.fs = gen_at;
      }
      d_load_check[1] = 1;
    }

    void parse_usb(Json::Value& signal){
      d_Sig_Params.am_norm = true;
      d_Sig_Params.mod_idx = signal.get("Modulation Index", 0.5).asFloat();
      d_Sig_Params.components = signal.get("Spectral Shaping Components", 0).asUInt();
      d_Sig_Params.spectral_len = signal.get("Spectral Shaping Tap Length", 256).asUInt();
      d_Sig_Params.f_max = signal.get("Spectral Shaping Max Frequency", 48e3).asFloat();
      d_Sig_Params.mu = std::vector<float>(d_Sig_Params.components,0.);
      d_Sig_Params.sigma = std::vector<float>(d_Sig_Params.components,0.);
      d_Sig_Params.weight = std::vector<float>(d_Sig_Params.components,0.);
      size_t idx = 0;
      for(Json::Value::iterator it1 = signal["Spectral Shaping Means"].begin(),
          it2 = signal["Spectral Shaping Std"].begin(),
          it3 = signal["Spectral Shaping Weights"].begin();
          ((it1 != signal["Spectral Shaping Means"].end()) &&
          (it2 != signal["Spectral Shaping Std"].end()) &&
          (it3 != signal["Spectral Shaping Weights"].end())) &&
          (idx < d_Sig_Params.components);
          it1++, it2++, it3++, idx++){
        d_Sig_Params.mu[idx] = (*it1).asFloat();
        d_Sig_Params.sigma[idx] = (*it2).asFloat();
        d_Sig_Params.weight[idx] = (*it3).asFloat();
      }
      d_Sig_Params.sps = signal.get("Interpolation", 1).asInt();
      d_Sig_Params.pulse_len = signal.get("Interpolation Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      idx = 0;
      for(Json::Value::iterator it = signal["Interpolation Taps"].begin();
          (it != signal["Interpolation Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      float gen_at = d_Sig_Params.f_max*d_Sig_Params.sps;
      float sanity = (d_Sig_Params.fs - gen_at)/gen_at;
      if(!((sanity > -1.19e-7)&&(sanity < 1.19e-7))){
        std::cerr << "signal_exciter:: json_parser:: signal_parameters:: " <<
            "parse_usb:: (WARN) 'Generation Sample Rate' != 'parameter space' " <<
            "defaulting to 'Generation Sample Rate' <= 'parameter space'" << std::endl;
        d_Sig_Params.fs = gen_at;
      }
      d_load_check[1] = 1;
    }

    void parse_lsb(Json::Value& signal){
      d_Sig_Params.am_norm = true;
      d_Sig_Params.mod_idx = signal.get("Modulation Index", 0.5).asFloat();
      d_Sig_Params.components = signal.get("Spectral Shaping Components", 0).asUInt();
      d_Sig_Params.spectral_len = signal.get("Spectral Shaping Tap Length", 256).asUInt();
      d_Sig_Params.f_max = signal.get("Spectral Shaping Max Frequency", 48e3).asFloat();
      d_Sig_Params.mu = std::vector<float>(d_Sig_Params.components,0.);
      d_Sig_Params.sigma = std::vector<float>(d_Sig_Params.components,0.);
      d_Sig_Params.weight = std::vector<float>(d_Sig_Params.components,0.);
      size_t idx = 0;
      for(Json::Value::iterator it1 = signal["Spectral Shaping Means"].begin(),
          it2 = signal["Spectral Shaping Std"].begin(),
          it3 = signal["Spectral Shaping Weights"].begin();
          ((it1 != signal["Spectral Shaping Means"].end()) &&
          (it2 != signal["Spectral Shaping Std"].end()) &&
          (it3 != signal["Spectral Shaping Weights"].end())) &&
          (idx < d_Sig_Params.components);
          it1++, it2++, it3++, idx++){
        d_Sig_Params.mu[idx] = (*it1).asFloat();
        d_Sig_Params.sigma[idx] = (*it2).asFloat();
        d_Sig_Params.weight[idx] = (*it3).asFloat();
      }
      d_Sig_Params.sps = signal.get("Interpolation", 1).asInt();
      d_Sig_Params.pulse_len = signal.get("Interpolation Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      idx = 0;
      for(Json::Value::iterator it = signal["Interpolation Taps"].begin();
          (it != signal["Interpolation Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      float gen_at = d_Sig_Params.f_max*d_Sig_Params.sps;
      float sanity = (d_Sig_Params.fs - gen_at)/gen_at;
      if(!((sanity > -1.19e-7)&&(sanity < 1.19e-7))){
        std::cerr << "signal_exciter:: json_parser:: signal_parameters:: " <<
            "parse_lsb:: (WARN) 'Generation Sample Rate' != 'parameter space' " <<
            "defaulting to 'Generation Sample Rate' <= 'parameter space'" << std::endl;
        d_Sig_Params.fs = gen_at;
      }
      d_load_check[1] = 1;
    }

    void parse_fm(Json::Value& signal){
      d_Sig_Params.mod_idx = signal.get("Modulation Index", 0.5).asFloat();
      d_Sig_Params.components = signal.get("Spectral Shaping Components", 0).asUInt();
      d_Sig_Params.spectral_len = signal.get("Spectral Shaping Tap Length", 256).asUInt();
      d_Sig_Params.f_max = signal.get("Spectral Shaping Max Frequency", 48e3).asFloat();
      d_Sig_Params.mu = std::vector<float>(d_Sig_Params.components,0.);
      d_Sig_Params.sigma = std::vector<float>(d_Sig_Params.components,0.);
      d_Sig_Params.weight = std::vector<float>(d_Sig_Params.components,0.);
      size_t idx = 0;
      for(Json::Value::iterator it1 = signal["Spectral Shaping Means"].begin(),
          it2 = signal["Spectral Shaping Std"].begin(),
          it3 = signal["Spectral Shaping Weights"].begin();
          ((it1 != signal["Spectral Shaping Means"].end()) &&
          (it2 != signal["Spectral Shaping Std"].end()) &&
          (it3 != signal["Spectral Shaping Weights"].end())) &&
          (idx < d_Sig_Params.components);
          it1++, it2++, it3++, idx++){
        d_Sig_Params.mu[idx] = (*it1).asFloat();
        d_Sig_Params.sigma[idx] = (*it2).asFloat();
        d_Sig_Params.weight[idx] = (*it3).asFloat();
      }
      d_Sig_Params.sps = signal.get("Interpolation", 1).asInt();
      d_Sig_Params.pulse_len = signal.get("Interpolation Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      idx = 0;
      for(Json::Value::iterator it = signal["Interpolation Taps"].begin();
          (it != signal["Interpolation Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      float gen_at = d_Sig_Params.f_max*d_Sig_Params.sps;
      float sanity = (d_Sig_Params.fs - gen_at)/gen_at;
      if(!((sanity > -1.19e-7)&&(sanity < 1.19e-7))){
        std::cerr << "signal_exciter:: json_parser:: signal_parameters:: " <<
            "parse_fm:: (WARN) 'Generation Sample Rate' != 'parameter space' " <<
            "defaulting to 'Generation Sample Rate' <= 'parameter space'" << std::endl;
        d_Sig_Params.fs = gen_at;
      }
      d_load_check[1] = 1;
      /*std::cout << "FM Loaded\n"
                << "mod idx " << d_Sig_Params.mod_idx  << std::endl
                << "comp " << d_Sig_Params.components  << std::endl
                << "spec len " << d_Sig_Params.spectral_len  << std::endl
                << "f max " << d_Sig_Params.f_max  << std::endl
                << "mu " << d_Sig_Params.mu.size()  << std::endl
                << "std " << d_Sig_Params.sigma.size()  << std::endl
                << "weight " << d_Sig_Params.weight.size()  << std::endl
                << "intp " << d_Sig_Params.sps  << std::endl
                << "intp len " << d_Sig_Params.pulse_len  << std::endl
                << "intp taps " << d_Sig_Params.pulse_shape.size()  << std::endl
                << std::endl;*/

    }

    void parse_psk(Json::Value& signal){
      d_Sig_Params.order = signal.get("Modulation Order", 2).asInt();
      d_Sig_Params.offset = signal.get("Phase Offset", 0.).asFloat();
      d_Sig_Params.sps = signal.get("Samples per Symbol", 2).asInt();
      d_Sig_Params.pulse_len = signal.get("Pulse Shape Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Pulse Shape Taps"].begin();
          (it != signal["Pulse Shape Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
    }

    void parse_qam(Json::Value& signal){
      d_Sig_Params.order = signal.get("Modulation Order", 2).asInt();
      d_Sig_Params.offset = signal.get("Phase Offset", 0.).asFloat();
      d_Sig_Params.sps = signal.get("Samples per Symbol", 2).asInt();
      d_Sig_Params.pulse_len = signal.get("Pulse Shape Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Pulse Shape Taps"].begin();
          (it != signal["Pulse Shape Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
    }

    void parse_pam(Json::Value& signal){
      d_Sig_Params.order = signal.get("Modulation Order", 2).asInt();
      d_Sig_Params.offset = signal.get("Phase Offset", 0.).asFloat();
      d_Sig_Params.sps = signal.get("Samples per Symbol", 2).asInt();
      d_Sig_Params.pulse_len = signal.get("Pulse Shape Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Pulse Shape Taps"].begin();
          (it != signal["Pulse Shape Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
    }

    void parse_ask(Json::Value& signal){
      d_Sig_Params.order = signal.get("Modulation Order", 2).asInt();
      d_Sig_Params.offset = signal.get("Phase Offset", 0.).asFloat();
      d_Sig_Params.sps = signal.get("Samples per Symbol", 2).asInt();
      d_Sig_Params.pulse_len = signal.get("Pulse Shape Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Pulse Shape Taps"].begin();
          (it != signal["Pulse Shape Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
    }

    void parse_ofdm(Json::Value& signal){
      std::string type = signal.get("Subcarrier Modulation", "").asString();
      std::transform( type.begin(), type.end(), type.begin(), ::tolower );
      if(strcmp(type.c_str(),"psk")==0){
        d_Sig_Params.mod = gr::signal_exciter::PSK;
      }
      else if(strcmp(type.c_str(), "qam")==0){
        d_Sig_Params.mod = gr::signal_exciter::QAM;
      }
      else if(strcmp(type.c_str(), "pam")==0){
        d_Sig_Params.mod = gr::signal_exciter::PAM;
      }
      else if(strcmp(type.c_str(), "ask")==0){
        d_Sig_Params.mod = gr::signal_exciter::ASK;
      }
      else{
        //unknown signal type
        std::cerr << "Signal_Exciter:: JSON_parser:: Signal_Parameters::" <<
          "\"Subcarrier Modulation\":: (WARN) Invalid Modulation Type -> " <<
          "Defaulting to PSK" << std::endl;
        d_Sig_Params.mod = gr::signal_exciter::PSK;
      }
      d_Sig_Params.order = signal.get("Modulation Order", 2).asInt();
      d_Sig_Params.offset = signal.get("Phase Offset", 0.).asFloat();
      d_Sig_Params.active_carriers = signal.get("Number of Subcarriers", 0).asUInt();
      d_Sig_Params.fftsize = signal.get("FFT Size", 0).asUInt();
      d_Sig_Params.cp_len = signal.get("Cyclic Prefix Length", 0).asUInt();
      d_Sig_Params.add_sync = signal.get("Add Schmidl and Cox Preamble", false).asBool();
      d_Sig_Params.syms_per_frame = signal.get("Symbols per Frame", 1).asUInt();
      d_Sig_Params.pilot_per_frame = signal.get("Pilots given per frame", false).asBool();
      d_Sig_Params.pilot_count = signal.get("Pilot Locations", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pilot_locations = std::vector<size_t>(d_Sig_Params.pilot_count,0);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Pilot Locations"].begin();
          (it != signal["Pilot Locations"].end()) &&
          (idx < d_Sig_Params.pilot_count);
          it++, idx++){
        d_Sig_Params.pilot_locations[idx] = (*it).asUInt();
      }
      d_Sig_Params.backoff = signal.get("Digital Backoff", 0.).asFloat();
      d_Sig_Params.samp_overlap = signal.get("Symbol Taper", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.taper = std::vector<float>(d_Sig_Params.samp_overlap);
      idx = 0;
      for(Json::Value::iterator it = signal["Symbol Taper"].begin();
          (it != signal["Symbol Taper"].end()) &&
          (idx < d_Sig_Params.samp_overlap);
          it++, idx++){
        d_Sig_Params.taper[idx] = (*it).asFloat();
      }
      d_Sig_Params.sps = signal.get("Interpolation", 1).asInt();
      d_Sig_Params.pulse_len = signal.get("Interpolation Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      idx = 0;
      for(Json::Value::iterator it = signal["Interpolation Taps"].begin();
          (it != signal["Interpolation Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
      std::cout << "OFDM Loaded\n"
                << "type " << type  << std::endl
                << "order " << d_Sig_Params.order  << std::endl
                << "offset " << d_Sig_Params.offset  << std::endl
                << "active " << d_Sig_Params.active_carriers  << std::endl
                << "fft " << d_Sig_Params.fftsize  << std::endl
                << "cp len " << d_Sig_Params.cp_len  << std::endl
                << "sync " << d_Sig_Params.add_sync  << std::endl
                << "ppf " << d_Sig_Params.pilot_per_frame  << std::endl
                << "pilot count " << d_Sig_Params.pilot_count  << std::endl
                << "pilots " << d_Sig_Params.pilot_locations.size()  << std::endl
                << "backoff " << d_Sig_Params.backoff << std::endl
                << "samp overlap " << d_Sig_Params.samp_overlap  << std::endl
                << "taper " << d_Sig_Params.taper.size()  << std::endl
                << "intp " << d_Sig_Params.sps  << std::endl
                << "intp len " << d_Sig_Params.pulse_len  << std::endl
                << "intp taps " << d_Sig_Params.pulse_shape.size()  << std::endl
                << std::endl;
    }

    void parse_ncfsk(Json::Value& signal){
      d_Sig_Params.order = signal.get("Modulation Order", 2).asInt();
      d_Sig_Params.sps = signal.get("Samples per Symbol", 2).asInt();
      d_Sig_Params.mod_idx = signal.get("Modulation Index", 0.5).asFloat();
      d_Sig_Params.pulse_len = signal.get("Augment Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Augment Taps"].begin();
          (it != signal["Augment Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
    }

    void parse_fsk(Json::Value& signal){
      d_Sig_Params.order = signal.get("Modulation Order", 2).asInt();
      d_Sig_Params.sps = signal.get("Samples per Symbol", 2).asInt();
      d_Sig_Params.mod_idx = signal.get("Modulation Index", 0.5).asFloat();
      d_Sig_Params.pulse_len = signal.get("Augment Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Augment Taps"].begin();
          (it != signal["Augment Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
    }

    void parse_cpm(Json::Value& signal){
      d_Sig_Params.order = signal.get("Modulation Order", 2).asInt();
      std::string type = signal.get("Phase Shape Type", "").asString();
      std::transform( type.begin(), type.end(), type.begin(), ::tolower );
      if(strcmp(type.c_str(),"lrec")==0){
        d_Sig_Params.phase_type = gr::analog::cpm::LREC;
      }
      else if(strcmp(type.c_str(), "lrc")==0){
        d_Sig_Params.phase_type = gr::analog::cpm::LRC;
      }
      else if(strcmp(type.c_str(), "lsrc")==0){
        d_Sig_Params.phase_type = gr::analog::cpm::LSRC;
      }
      else if(strcmp(type.c_str(), "tfm")==0){
        d_Sig_Params.phase_type = gr::analog::cpm::TFM;
      }
      else if(strcmp(type.c_str(), "gaussian")==0){
        d_Sig_Params.phase_type = gr::analog::cpm::GAUSSIAN;
      }
      else{
        //unknown signal type
        std::cerr << "Signal_Exciter:: JSON_parser:: Signal_Parameters::" <<
          "\"Phase Shape Type\":: (WARN) Invalid Phase Shape Type -> " <<
          "Defaulting to LREC" << std::endl;
        d_Sig_Params.phase_type = gr::analog::cpm::LREC;
      }
      d_Sig_Params.sps = signal.get("Samples per Symbol", 2).asInt();
      d_Sig_Params.L = signal.get("Symbol Overlap", 1).asInt();
      d_Sig_Params.mod_idx = signal.get("Modulation Index", 0.5).asFloat();
      d_Sig_Params.beta = signal.get("Phase Shape Beta", 0.35).asFloat();
      d_Sig_Params.pulse_len = signal.get("Augment Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Augment Taps"].begin();
          (it != signal["Augment Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
    }

    void parse_gmsk(Json::Value& signal){
      d_Sig_Params.sps = signal.get("Samples per Symbol", 2).asInt();
      d_Sig_Params.L = signal.get("Symbol Overlap", 1).asInt();
      d_Sig_Params.beta = signal.get("Phase Shape Beta", 0.35).asFloat();
      d_Sig_Params.pulse_len = signal.get("Augment Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Augment Taps"].begin();
          (it != signal["Augment Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
    }

    void parse_msk(Json::Value& signal){
      d_Sig_Params.sps = signal.get("Samples per Symbol", 2).asInt();
      d_Sig_Params.pulse_len = signal.get("Augment Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Augment Taps"].begin();
          (it != signal["Augment Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
    }

    void parse_gfsk(Json::Value& signal){
      d_Sig_Params.order = signal.get("Modulation Order", 2).asInt();
      d_Sig_Params.sps = signal.get("Samples per Symbol", 2).asInt();
      d_Sig_Params.L = signal.get("Symbol Overlap", 1).asInt();
      d_Sig_Params.mod_idx = signal.get("Modulation Index", 0.5).asFloat();
      d_Sig_Params.beta = signal.get("Phase Shape Beta", 0.35).asFloat();
      d_Sig_Params.pulse_len = signal.get("Augment Taps", Json::Value(Json::arrayValue)).size();
      d_Sig_Params.pulse_shape = std::vector<float>(d_Sig_Params.pulse_len,0.);
      size_t idx = 0;
      for(Json::Value::iterator it = signal["Augment Taps"].begin();
          (it != signal["Augment Taps"].end()) &&
          (idx < d_Sig_Params.pulse_len);
          it++, idx++){
        d_Sig_Params.pulse_shape[idx] = (*it).asFloat();
      }
      d_load_check[1] = 1;
    }

    void parse_gating(Json::Value& signal){
      d_Sig_Params.ops_gate = signal.get("One Pass Gate", false).asBool();
      d_Sig_Params.per_gate = signal.get("Periodic Gate", false).asBool();
      d_Sig_Params.rnd_gate = signal.get("Random Gate", false).asBool();
      d_Sig_Params.ops_gate_off = signal.get("One Pass Gate Off Duration", 0.).asFloat();
      d_Sig_Params.ops_gate_on = signal.get("One Pass Gate On Duration", 0.).asFloat();
      d_Sig_Params.per_gate_off = signal.get("Periodic Gate Off Duration", 0.).asFloat();
      d_Sig_Params.per_gate_on = signal.get("Periodic Gate On Duration", 0.).asFloat();
      d_Sig_Params.per_gate_offset = signal.get("Periodic Gate Period Offset", 0.).asFloat();
      d_Sig_Params.rnd_gate_off_min = signal.get("Random Gate Off Duration Min", 0.).asFloat();
      d_Sig_Params.rnd_gate_off_max = signal.get("Random Gate Off Duration Max", 0.).asFloat();
      d_Sig_Params.rnd_gate_on_min = signal.get("Random Gate Off Duration Min", 0.).asFloat();
      d_Sig_Params.rnd_gate_on_max = signal.get("Random Gate Off Duration Max", 0.).asFloat();
      d_load_check[2] = 1;
    }

    void parse_general(Json::Value& signal){
      bool loaded = true;
      d_Sig_Params.fc = signal.get("Center Frequency", 0.).asFloat();
      d_Sig_Params.fs = signal.get("Generation Sample Rate", 0.).asFloat();
      d_Sig_Params.gain = signal.get("Gain", 0.).asFloat();
      Json::Value frac_check = signal.get("Fractional Symbol Offset", Json::Value());
      if(frac_check.isNull()){
        d_Sig_Params.frac_offset = false;
        d_Sig_Params.frac_symb_offset = 0.;
      }
      else{
        d_Sig_Params.frac_offset = true;
        d_Sig_Params.frac_symb_offset = signal.get("Fractional Symbol Offset", 0.).asFloat();
      }
      d_load_check[0] = 1;
      std::cout << "General loaded\n" << d_Sig_Params.fc << std::endl
                << d_Sig_Params.fs << std::endl
                << d_Sig_Params.gain << std::endl
                << d_Sig_Params.frac_symb_offset << std::endl;
      std::string type = signal.get("Modulation Type", "").asString();
      std::transform( type.begin(), type.end(), type.begin(), ::tolower );
      if(strcmp(type.c_str(),"cwmorse")==0){
        std::cout << "CWMORSE\n";
        d_Sig_Params.type = gr::signal_exciter::CWMORSE;
        parse_cwmorse(signal);
      }
      else if(strcmp(type.c_str(), "dsb")==0){
        std::cout << "DSB\n";
        d_Sig_Params.type = gr::signal_exciter::DSB;
        parse_dsb(signal);
      }
      else if(strcmp(type.c_str(), "dsbsc")==0){
        std::cout << "DSBSC\n";
        d_Sig_Params.type = gr::signal_exciter::DSBSC;
        parse_dsbsc(signal);
      }
      else if(strcmp(type.c_str(), "usb")==0){
        std::cout << "USB\n";
        d_Sig_Params.type = gr::signal_exciter::USB;
        parse_usb(signal);
      }
      else if(strcmp(type.c_str(), "lsb")==0){
        std::cout << "LSB\n";
        d_Sig_Params.type = gr::signal_exciter::LSB;
        parse_lsb(signal);
      }
      else if(strcmp(type.c_str(), "fm")==0){
        std::cout << "FM\n";
        d_Sig_Params.type = gr::signal_exciter::FM;
        std::cout << "Initializing FM Signal" << std::endl;
        parse_fm(signal);
      }
      else if(strcmp(type.c_str(), "psk")==0){
        std::cout << "PSK\n";
        d_Sig_Params.type = gr::signal_exciter::PSK;
        parse_psk(signal);
      }
      else if(strcmp(type.c_str(), "qam")==0){
        std::cout << "QAM\n";
        d_Sig_Params.type = gr::signal_exciter::QAM;
        parse_qam(signal);
      }
      else if(strcmp(type.c_str(), "pam")==0){
        std::cout << "PAM\n";
        d_Sig_Params.type = gr::signal_exciter::PAM;
        parse_pam(signal);
      }
      else if(strcmp(type.c_str(), "ask")==0){
        std::cout << "ASK\n";
        d_Sig_Params.type = gr::signal_exciter::ASK;
        parse_ask(signal);
      }
      else if(strcmp(type.c_str(), "ofdm")==0){
        std::cout << "OFDM\n";
        d_Sig_Params.type = gr::signal_exciter::OFDM;
        parse_ofdm(signal);
      }
      else if(strcmp(type.c_str(), "ncfsk")==0){
        std::cout << "NCFSK\n";
        d_Sig_Params.type = gr::signal_exciter::NCFSK;
        parse_fsk(signal);
      }
      else if(strcmp(type.c_str(), "fsk")==0){
        std::cout << "FSK\n";
        d_Sig_Params.type = gr::signal_exciter::FSK;
        parse_fsk(signal);
      }
      else if(strcmp(type.c_str(), "cpm")==0){
        std::cout << "CPM\n";
        d_Sig_Params.type = gr::signal_exciter::CPM;
        parse_cpm(signal);
      }
      else if(strcmp(type.c_str(), "gmsk")==0){
        std::cout << "GMSK\n";
        d_Sig_Params.type = gr::signal_exciter::GMSK;
        parse_gmsk(signal);
      }
      else if(strcmp(type.c_str(), "msk")==0){
        std::cout << "MSK\n";
        d_Sig_Params.type = gr::signal_exciter::MSK;
        parse_msk(signal);
      }
      else if(strcmp(type.c_str(), "gfsk")==0){
        std::cout << "GFSK\n";
        d_Sig_Params.type = gr::signal_exciter::GFSK;
        parse_gfsk(signal);
      }
      else{
        //unknown signal type
        std::cerr << "Signal_Exciter:: JSON_parser:: Signal_Parameters::" <<
          "\"Modulation Type\":: (WARN) Unknown Modulation Type (skip)" << std::endl;
      }
      load_check();
    }

   public:
    /*
     *  \brief Empty Signal_Parameters Initialize
     *
     *
     */
    Signal_Parameters(){
      d_Loaded = false;
      d_load_check = std::vector<char>(3,0);
    }

    /*
     *  \brief Initalize Parameters from a JSON parser object
     *
     *
     */
    Signal_Parameters(Json::Value& signal)
    {
      std::cout << "INIT SigParam\n";
      d_Loaded = false;
      d_load_check = std::vector<char>(3,0);
      parse_general(signal);
      parse_gating(signal);
    }

    bool is_loaded(){
      if(!d_Loaded){
        load_check();
      }
      return d_Loaded;
    }
    bool get_params(gr::signal_exciter::sig_params &sig_param){
      if(d_Loaded){
        sig_param = d_Sig_Params;
        return true;
      }
      return false;
    }

    Json::Value get_json(){
      if(d_Loaded){
        return Json::Value(Json::objectValue);
      }
    }

  };

  class System_Parameters
  {
   private:
    int               d_Global_Seed;
    size_t            d_Signal_Count;
    double            d_System_Center;
    double            d_System_Bandwidth;
    double            d_Run_Time;
    double            d_Noise_Floor;
    bool              d_Noise_Enabled;
    std::string       d_File_Location;
    bool              d_Loaded;

    std::vector<Signal_Parameters> d_Signals;

    std::vector<char> d_load_check;
    void load_check(){
      bool loaded = true;
      size_t idx = 0;
      while(loaded && (idx < 9)){
        loaded = loaded && d_load_check[idx++];
      }
      if(loaded && (idx>=9)){ d_Loaded = true; }
    }

    std::vector<char> d_signal_check;
    bool signal_check(){
      bool loaded = true;
      size_t idx = 0;
      while(loaded && (idx < d_Signals.size())){
        loaded = loaded && d_signal_check[idx++];
      }
      return loaded;
    }

   public:
    /*
     *  \brief Empty System_Parameters Initialize
     *
     *
     */
    System_Parameters()
    {
      d_Loaded = false;
      d_Signal_Count = 0;
      d_load_check = std::vector<char>(9,0.);
    }

    /*
     *  \brief Initialize System_Parameters with File_Location
     *
     *
     */
    System_Parameters(const std::string& File_Location)
    {
      d_File_Location = std::string(File_Location);
      d_Loaded = false;
      d_Signal_Count = 0;
      d_load_check = std::vector<char>(9,0.);
      d_load_check[0] = 1;
    }

    ~System_Parameters(){}

    int get_next_seed()
    {
      if(d_Global_Seed < 0){
        return d_Global_Seed;
      }
      d_Global_Seed++;
      if(d_Global_Seed < 0) d_Global_Seed = 0;
      return d_Global_Seed;
    }

    bool is_loaded(){
      if(!d_Loaded){
        signal_check();
        load_check();
      }
      return d_Loaded;
    }
    int get_seed(){ return d_Global_Seed; }
    size_t get_signal_count(){ return d_Signal_Count; }
    double get_system_center(){ return d_System_Center; }
    double get_system_bandwidth(){ return d_System_Bandwidth; }
    double get_run_time(){ return d_Run_Time; }
    double get_noise_floor(){ return d_Noise_Floor; }
    bool get_noise_enabled(){ return d_Noise_Enabled; }
    std::string get_file_location(){ return d_File_Location; }
    std::vector<Signal_Parameters>& get_signal_parameters(){ return d_Signals; }
    bool get_signal_parameters(gr::signal_exciter::sig_params &sp, size_t idx){
      if(idx < d_Signals.size()){
        return d_Signals[idx].get_params(sp);
      }
      return false;
    }

    void set_seed(int seed)
    {
      d_Global_Seed = seed;
      d_load_check[1] = 1;
      is_loaded();
      std::cout<<"Loaded Seed\n";
    }
    void set_signal_count(size_t signal_count){
      d_Signal_Count = signal_count;
      d_Signals = std::vector<Signal_Parameters>(signal_count);
      d_signal_check = std::vector<char>(signal_count,0);
      d_load_check[2] = 1;
      is_loaded();
      std::cout<<"Loaded SC " << d_Signal_Count << "\n";
    }
    void set_system_center(double sys_fc){
      d_System_Center = sys_fc;
      d_load_check[3] = 1;
      is_loaded();
      std::cout<<"Loaded fc " << d_System_Center << "\n";
    }
    void set_system_bandwidth(double sys_bw){
      d_System_Bandwidth = sys_bw;
      d_load_check[4] = 1;
      is_loaded();
      std::cout<<"Loaded fs " << d_System_Bandwidth << "\n";
    }
    void set_run_time(double run_time){
      d_Run_Time = run_time;
      d_load_check[5] = 1;
      is_loaded();
      std::cout<<"Loaded rt " << d_Run_Time << "\n";
    }
    void set_noise_enable(bool enable){
      d_Noise_Enabled = enable;
      d_load_check[6] = 1;
      is_loaded();
      std::cout<<"Loaded ne " << d_Noise_Enabled << "\n";
    }
    void set_noise_floor(double noise_floor){
      d_Noise_Floor = noise_floor;
      d_load_check[7] = 1;
      is_loaded();
      std::cout<<"Loaded nf " << d_Noise_Floor << "\n";
    }
    void set_file_location(const std::string& File_Location){
      d_File_Location = std::string(File_Location);
      d_load_check[0] = 1;
      is_loaded();
      std::cout<<"Loaded fl " << d_File_Location << "\n";
    }
    void set_signal_parameters(std::vector<Signal_Parameters>& Signals){
      d_Signals = std::vector<Signal_Parameters>(Signals.begin(), Signals.end());
      if(d_Signals.size() == d_Signal_Count){
        d_load_check[8] = 1;
        is_loaded();
      }
    }
    bool set_signal_parameters(Signal_Parameters& Signal, size_t idx){
      if(d_load_check[2]){
        if(idx < d_Signal_Count){
          d_Signals[idx] = Signal;
          d_signal_check[idx] = 1;
          is_loaded();
          std::cout<<"Loaded signal " << idx << std::endl ;
        }
        else{
          throw std::out_of_range("Signal_Exciter:: JSON_parser:: System_Parameters::"
            "set_signal_parameters:: index location out of range\n");
        }
      }
      else{
        throw std::runtime_error("Signal_Exciter:: JSON_parser:: System_Parameters::"
          "set_signal_parameters:: set_signal_count must be called first\n");
      }
      d_load_check[8] = signal_check();
      return true;
    }

    /*
     *  \brief Read in the provided file name
     *
     *
     */
    bool Load_File(){
      if(d_load_check[0]){
        std::ifstream infile(d_File_Location.c_str(), std::ifstream::binary);
        Json::Value json_book;
        infile >> json_book;
        set_system_bandwidth(json_book.get("System Bandwidth", 1.).asDouble());
        set_system_center(json_book.get("Center Frequency", 0.).asDouble());
        set_signal_count(json_book.get("Number of Signals", 0).asUInt());
        set_run_time(json_book.get("Run Time", 0.).asDouble());
        set_noise_enable(json_book.get("Noise Enabled", false).asBool());
        set_noise_floor(json_book.get("Noise Floor", -100.).asDouble());
        set_file_location(json_book.get("Output File", "something.fc32").asString());
        set_seed(json_book.get("Global Seed", -1).asInt());

        bool good_load = true;
        size_t idx = 0;
        for(Json::Value::iterator it = json_book["Signals"].begin();
            it != json_book["Signals"].end();
            it++){
          idx = (*it).get("ID", int(get_signal_count())).asUInt();
          std::cout << "Testing " << idx << std::endl;
          if(idx < get_signal_count()){
            Signal_Parameters sig((*it));
            if(sig.is_loaded()){
              set_signal_parameters(sig, idx);
            }
            else{
              good_load = false;
              std::cerr << "Signal " << idx << " failed to load" << std::endl;
            }
          }
        }
        return good_load;
      }
      return false;
    }


    /*
     *  \brief Save the parameters to the provided file name
     *
     *
     */
    bool Save_File(){

    }

  };

}//json
}//signal_exciter

#endif /* INCLUDED_SIGNAL_EXCITER_JSON_PARSER_HPP */
