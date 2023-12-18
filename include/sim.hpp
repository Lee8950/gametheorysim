#ifndef ECL_SIM_HPP
#define ECL_SIM_HPP 1

#include <bits/stdc++.h>

namespace ecl
{

class zipf {
public:
    double sum;

    zipf(int K) {
        sum = 0;
        for(int i = 1; i <= K+1; i++)
            sum += 1.0/i;
    }

    double operator()(int k) {
        return 1.0 / (sum * k);
    }

};

extern int S;
extern int Q;

extern double c0;
extern double cs;
extern double epsilon;

extern int file_count;
extern int device_count;

extern std::random_device rd;
extern std::mt19937 gen;
extern std::uniform_real_distribution<> uniform;
extern std::vector<double> tracker;

class Device;
class File;

class Device {

public:

    double total_demand;
    double caching_cost;
    double serving_capacity;
    double storage_capacity;

    int device_id;

public:

    std::vector<int> caching_desicions;

    Device(double total_demand, double caching_cost, double serving_capacity, double storage_capacity, int device_id, int files_total):
        total_demand(total_demand), caching_cost(caching_cost), serving_capacity(serving_capacity), storage_capacity(storage_capacity), device_id(device_id) {
        
        caching_desicions.resize(files_total);
        memset(caching_desicions.data(), 0, sizeof(int) * files_total);
    }

    

public:

    double U1(int k, std::vector<ecl::File>& files);

    double U0(int k, std::vector<ecl::File>& files);

    double Us(int k, std::vector<ecl::File>& files);

    double beta(int k, std::vector<ecl::File>& files);


public:

    // Algorithm 1
    std::vector<int> MEC_PEC_Strategy_Update(std::vector<double> alphas, std::vector<double> fis, std::vector<File>& files, std::vector<Device>& devices);

};

class File {

private:

public:

    double v;
    double popularity;
    double price;

public:

    File(double v, double popularity, double price):
        v(v), popularity(popularity), price(price) {}

};

std::vector<double> pi1(int device_index, std::vector<ecl::File>& files, std::vector<ecl::Device>& devices, std::vector<double>& fis);

std::vector<double> pi0(int device_index, std::vector<ecl::File>& files, std::vector<ecl::Device>& devices, std::vector<double>& alphas);

double Dk(int k, std::vector<File>& files, std::vector<Device>& devices);

int S1Count(int k, std::vector<Device>& devices);

int S0Count(int k, std::vector<Device>& devices);

int iterativeS1(int k, std::vector<File>& files, std::vector<Device>& devices, std::vector<double>& alphas, std::vector<double>& fis);

int iterativeS0(int k, std::vector<File>& files, std::vector<Device>& devices, std::vector<double>& alphas, std::vector<double>& fis);

double get_alpha(int k, double Q, std::vector<File>& files, std::vector<Device>& devices);

double get_fi(int k, double Q, std::vector<File>& files, std::vector<Device>& devices);

}

#endif