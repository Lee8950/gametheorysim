#include <sim.hpp>
#include <nlohmann/json.hpp>

int main(int argc, char** argv)
{
    std::string prefix;
    
    // Set hyper params
    ecl::device_count = 2000;
    ecl::file_count = 1000;
    ecl::epsilon = 1e-8;
    ecl::cs = 0.4e-3;
    ecl::Q = 1000;
    ecl::S = 20;
    ecl::c0 = 0.2e-3;
    
    double demand_factor = 1;
    double discount_factor = 0;
    double constant_sharing_price = 0.1e-3;
    double update_rate = 0.1;
    if(argc>2)
    {
        prefix = argv[1];
        ecl::Q = atof(argv[2]);
        constant_sharing_price = atof(argv[3]);
        if(argc == 5)
            ecl::c0 = atof(argv[4]);
    }
    //std::cout << prefix << ecl::Q << constant_sharing_price;
    //return 0;
    ecl::zipf zipf_dist(ecl::file_count);
    double lambda = 0.2;
    ecl::debug_var = 0;

    // Create devices and files
    std::vector<ecl::Device> devices;
    std::vector<ecl::File> files;
    
    for(int i = 0; i < ecl::device_count; i++)
        devices.push_back(ecl::Device(demand_factor * ecl::wdist(ecl::gen), ecl::uniform(ecl::gen), ecl::Q, ecl::S, i, ecl::file_count));

    for(int i = 0; i < ecl::file_count; i++)
        files.push_back(ecl::File(0, zipf_dist(i+1), constant_sharing_price+discount_factor*zipf_dist(i+1)));


    // Compute init fi and alpha
    
    std::vector<double> fis;
    std::vector<double> alphas;
    fis.resize(ecl::file_count);
    alphas.resize(ecl::file_count);
    memset(fis.data(), 0, sizeof(double) * ecl::file_count);
    memset(alphas.data(), 0, sizeof(double) * ecl::file_count);
    {
        auto newfis = fis;
        auto newalphas = alphas;
        for(int i = 0; i < ecl::file_count; i++)
        {
            newalphas[i] = ecl::get_alpha(i, ecl::Q, files, devices);
            newfis[i] = lambda * fis[i] + (1 - lambda) * ecl::Dk(i, files, devices) * files[i].price / ecl::iterativeS1(i, files, devices, alphas, fis);
            //newfis[i] = ecl::get_fi(i, ecl::Q, files, devices);
            //newfis[i] = lambda * fis[i] + (1 - lambda) * ecl::Dk(i, files, devices) * files[i].price / ecl::S1Count(i, devices);
        }
        for(int i = 0; i < ecl::file_count; i++)
        {
            alphas[i] = newalphas[i];
            fis[i] = newfis[i];
        }
    }


    // while loop
    int test = 0;
    printf("\n");
    auto bar = [](int a, int b)->std::string{
        int barLength = 50;
        int stop = barLength * (static_cast<double>(a) / static_cast<double>(b));
        std::string s;
        for(int i = 0; i < stop; i++) s.append("=");
        for(int i = 0; i < barLength-stop; i++) s.append(".");
        return s;
    };
    
    // logging
    std::fstream loggingfile;
    std::fstream latest;
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::string logLocation; //std::put_time(std::localtime(&now), "%Y.%m.%d_%H:%M:%S");
    std::stringstream tts;
    tts << std::put_time(std::localtime(&now), "%Y.%m.%d_%H.%M.%S");
    if(argc == 1)
        logLocation = std::filesystem::absolute(argv[0]).parent_path() / "logs" / (tts.str() + ".log"); //"logs/" + tts.str() + ".log";
    else
        logLocation = std::filesystem::absolute(argv[0]).parent_path() / "logs" / (prefix + tts.str() + ".log");
    loggingfile.open(logLocation, std::ios::out);
    latest.open(std::filesystem::absolute(argv[0]).parent_path() / "logs/latest.log", std::ios::out);
    if(!loggingfile.is_open() || !latest.is_open())
        exit(EXIT_FAILURE);

    std::vector<int> f1s, f2s, f3s, f4s;

    int epoch = 100;
    int f1=0, f2=0, f3=0, f4=0;
    while(test++<epoch)
    {
        ecl::debug_var = test;
        printf("\033[1AEpoch:%d/%d <%s>\n",test, epoch, bar(test, epoch).c_str());
        std::vector<std::vector<int>> strategies;
        strategies.resize(ecl::device_count);
        for(int i = 0; i < ecl::device_count; i++)
            strategies[i] = devices[i].MEC_PEC_Strategy_Update(alphas, fis, files, devices);

        auto newfis = fis;
        auto newalphas = alphas;
        for(int i = 0; i < ecl::file_count; i++)
        {
            newalphas[i] = ecl::get_alpha(i, ecl::Q, files, devices);
            newfis[i] = lambda * fis[i] + (1 - lambda) * ecl::Dk(i, files, devices) * files[i].price / ecl::iterativeS1(i, files, devices, alphas, fis);
            //newfis[i] = ecl::get_fi(i, ecl::Q, files, devices);
            //newfis[i] = lambda * fis[i] + (1 - lambda) * ecl::Dk(i, files, devices) * files[i].price / ecl::S1Count(i, devices);
        }
        for(int i = 0; i < ecl::file_count; i++)
        {
            alphas[i] = update_rate * newalphas[i] + (1 - update_rate) * alphas[i];
            fis[i] = update_rate * newfis[i] + (1 - update_rate) * fis[i];
        }
        for(int i = 0; i < ecl::device_count; i++)
            devices[i].caching_desicions = strategies[i];
        // statistics
        f1=0, f2=0, f3=0, f4=0;
        for(int i = 0; i < ecl::device_count; i++)
        {
            f1 += (devices[i].caching_desicions[0] == 1);
            f2 += (devices[i].caching_desicions[1] == 1);
            f3 += (devices[i].caching_desicions[2] == 1);
            f4 += (devices[i].caching_desicions[3] == 1);
        }
        f1s.push_back(f1);
        f2s.push_back(f2);
        f3s.push_back(f3);
        f4s.push_back(f4);


    }

    f1=0, f2=0, f3=0, f4=0;
    for(int i = 0; i < ecl::device_count; i++)
    {
        f1 += (devices[i].caching_desicions[0] == 1);
        f2 += (devices[i].caching_desicions[1] == 1);
        f3 += (devices[i].caching_desicions[2] == 1);
        f4 += (devices[i].caching_desicions[3] == 1);
    }

    for(int i = 0; i < 10; i++)
    {
        int stat = 0;
        for(int j = 0; j < ecl::device_count; j++)
            stat += (devices[j].caching_desicions[i] == 1);
        printf("File %d Agent:%d/%d, %.2lf%%\n", i+1, stat, ecl::device_count, 100 * stat / static_cast<double>(ecl::device_count));
    }

    double system_cost = 0.0;
    for(int i = 0; i < ecl::device_count; i++)
    {
        for(int j = 0; j < ecl::file_count; j++)
        {
            if(devices[i].caching_desicions[j] == 1)
            {
                system_cost += devices[i].caching_cost;
            }
            else
            {
                double u0 = devices[i].U0(j, files);
                double beta = devices[i].beta(j, files);
                double us = devices[i].Us(j, files);
                system_cost += devices[i].total_demand * files[j].popularity * (alphas[j] * ecl::c0 + (1 - alphas[j]) * ecl::cs);
                //system_cost += alphas[j] * (u0 - beta) + (1-alphas[j])*us;
            }
        }
    }
    printf("system_cost = %.2lf\n", system_cost);

    /*printf("File 1 Agent:%d/%d, %.2lf%%\nFile 2 Agent:%d/%d, %.2lf%%\nFile 3 Agent:%d/%d, %.2lf%%\nFile 4 Agent:%d/%d, %.2lf%%\n",
        f1, ecl::device_count, 100 * f1 / static_cast<double>(ecl::device_count),
        f2, ecl::device_count, 100 * f2 / static_cast<double>(ecl::device_count),
        f3, ecl::device_count, 100 * f3 / static_cast<double>(ecl::device_count),
        f4, ecl::device_count, 100 * f4 / static_cast<double>(ecl::device_count));*/

    nlohmann::json target;
    target["device_count"] = ecl::device_count;
    target["file_count"] = ecl::file_count;
    target["cs"] = ecl::cs;
    target["Q"] = ecl::Q;
    target["S"] = ecl::S;
    target["c0"] = ecl::c0;
    target["demand_factor"] = demand_factor;
    target["discount_factor"] = discount_factor;
    target["constant_sharing_price"] = constant_sharing_price;
    target["epoch"] = epoch;
    target["f1"] = f1s;
    target["f2"] = f2s;
    target["f3"] = f3s;
    target["f4"] = f4s;
    target["system_cost_sum_files_and_devices"] = -system_cost;
    loggingfile << target.dump();
    latest << target.dump();
    latest.close();
    loggingfile.close();
    auto t = ecl::tracker;
    return 0;
}