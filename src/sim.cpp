#include <sim.hpp>

namespace ecl {

int S;
int Q;

double c0;
double cs;
double epsilon;

int file_count;
int device_count;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> uniform(0.0, 1.0);
std::vector<double> tracker;

}

double ecl::Device::U1(int k, std::vector<ecl::File> &files)
{
    return total_demand * files[k].popularity * files[k].v - caching_cost;
}

double ecl::Device::U0(int k, std::vector<ecl::File> &files)
{
    return total_demand * files[k].popularity * (files[k].v - c0);
}

double ecl::Device::Us(int k, std::vector<ecl::File> &files)
{
    return total_demand * files[k].popularity * (files[k].v - cs);
}

double ecl::Device::beta(int k, std::vector<ecl::File> &files)
{
    return total_demand * files[k].popularity * files[k].price;
}

double ecl::Dk(int k, std::vector<ecl::File> &files, std::vector<ecl::Device> &devices)
{
    double s = 0;
    double factor = 0;
    for(Device& device : devices)
    {
        s = s + files[k].popularity * device.total_demand;
        for(int selection : device.caching_desicions)
            factor += (selection == 0);
    }
    tracker.push_back(s/factor);
    return s / factor;
}

int ecl::S1Count(int k, std::vector<ecl::Device> &devices)
{
    int s = 0;
    for(Device& device : devices)
        s += (device.caching_desicions[k] == 1);
    return s;
}

int ecl::S0Count(int k, std::vector<ecl::Device> &devices)
{
    int s = 0;
    for(Device& device : devices)
        s += (device.caching_desicions[k] == 0);
    return s;
}

double ecl::get_alpha(int k, double Q, std::vector<ecl::File> &files, std::vector<ecl::Device> &devices)
{
    return std::min(S0Count(k, devices) * Q / Dk(k, files, devices), static_cast<double>(1.0));
}

double ecl::get_fi(int k, double Q, std::vector<ecl::File> &files, std::vector<ecl::Device> &devices)
{
    int s1 = S1Count(k, devices);
    if(s1 == 0)
        return 1;
    return get_alpha(k, Q, files, devices)*Dk(k, files, devices)*files[k].price/s1;
}

std::vector<double> ecl::pi1(int device_index, std::vector<ecl::File> &files, std::vector<ecl::Device> &devices, std::vector<double> &fis)
{
    std::vector<double> tmp;
    tmp.resize(files.size());
    for(int i = 0; i < tmp.size(); i++)
    {
        tmp[i] = devices[device_index].U1(i, files) + fis[i];
    }
    return tmp;
}

std::vector<double> ecl::pi0(int device_index, std::vector<ecl::File> &files, std::vector<ecl::Device> &devices, std::vector<double> &alphas)
{
    std::vector<double> tmp;
    tmp.resize(files.size());
    for(int i = 0; i < tmp.size(); i++)
    {
        tmp[i] = alphas[i] * (devices[device_index].U0(i, files) - devices[device_index].beta(i, files)) + 
        (1 - alphas[i]) * devices[device_index].Us(i, files);
    }
    return tmp;
}

std::vector<int> ecl::Device::MEC_PEC_Strategy_Update(std::vector<double> alphas, std::vector<double> fis, std::vector<ecl::File> &files, std::vector<ecl::Device> &devices)
{
    std::vector<int> new_caching_desicions;
        new_caching_desicions.resize(file_count);
        memset(new_caching_desicions.data(), 0, sizeof(int) * file_count);

        std::vector<double> devicePi0 = pi0(device_id, files, devices, alphas);
        std::vector<double> devicePi1 = pi1(device_id, files, devices, fis);
        std::vector<int> seq;
        seq.resize(devicePi0.size());
        for(int i = 0; i < seq.size(); i++) seq[i] = i;

        for(int i = 0; i < seq.size(); i++)
        {
            for(int j = 0; j < seq.size() - 1; j++)
            {
                if((devicePi1[j]-devicePi0[j]) < devicePi1[j+1]-devicePi0[j+1])
                {
                    std::swap(devicePi1[j], devicePi1[j+1]);
                    std::swap(devicePi0[j], devicePi0[j+1]);
                    std::swap(seq[j], seq[j+1]);
                }
            }
        }

        auto currentCapacity = [](std::vector<int> desicions)->int{
            int cap = 0;
            for(auto& selection : desicions)
            {
                cap += (selection == 1);
            }
            return cap;
        };

        for(int seq_index = 0; seq_index < files.size(); seq_index++)
        {
            if(devicePi1[seq[seq_index]] - devicePi0[seq[seq_index]] >= epsilon && currentCapacity(new_caching_desicions) + 1 <= storage_capacity)
                new_caching_desicions[seq[seq_index]] = 1;
            else
                new_caching_desicions[seq[seq_index]] = 0;
        }

        return new_caching_desicions;
}

int ecl::iterativeS1(int k, std::vector<File>& files, std::vector<ecl::Device> &devices, std::vector<double> &alphas, std::vector<double> &fis)
{
    int s = 0;
    for(Device& device : devices)
        if(device.caching_cost <= fis[k]+device.total_demand*files[k].popularity*(alphas[k]*(files[k].price+c0)+(1-alphas[k]*cs)))
            s++;
    return s;
}

int ecl::iterativeS0(int k, std::vector<File>& files, std::vector<ecl::Device> &devices, std::vector<double> &alphas, std::vector<double> &fis)
{
    int s = 0;
    for(Device& device : devices)
        if(device.caching_cost > fis[k]+device.total_demand*files[k].popularity*(alphas[k]*(files[k].price+c0)+(1-alphas[k]*cs)))
            s++;
    return s;
}