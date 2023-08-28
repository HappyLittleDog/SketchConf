#include "experiment.h"
using namespace std;

count_t TINY_THRESHOLD;
double BIG_PERCENT;

extern FILE *yyin;
extern int yyparse(std::unique_ptr<BaseAst> &ast);

string get_sketch_name(int t)
{
    string s;
    switch (t)
    {
    case int(SketchType::CM):
        s="CM";
        break;

    case int(SketchType::COUNT):
        s="Count";
        break;
    
    case int(SketchType::TOWER):
        s="Tower";
        break;
    
    case int(SketchType::NITROCM):
        s="Nitro CM";
        break;
    
    case int(SketchType::NITROCOUNT):
        s="Nitro Count";
        break;
    
    default:
        LOG_ERROR("Unexpected sketch type: %d", int(t));
        exit(-1);
    }
    return s;
}

void experiment_1()
{
    LOG_INFO("Experiment #1: Experiments on the Error Predictor");
    Constraint constraint;
    constraint.err=100;
    constraint.prob=1.0;
    deque<Constraint> constraints;
    constraints.push_back(constraint);
    StreamGen stream;
    stream.init("dataset/data.dat", 4);
    int fd;
    double sum;
    double round=10;

    LOG_INFO("Current Sketch: CM");
    fd=Open("output/exp1-cm.txt", O_WRONLY | O_CREAT);
    CountMinSketch cm(1);
    LOG_RESULT("Ground Truth: %lf", cm.ground_truth(3, 100000, constraints, stream)[0]);
    sum=0;
    for (int i=0;i<round;i++)
    {
        double ans=cm.simulate(3, 100000, constraints, stream, false, false, false)[0];
        sum += ans;
        string cur=to_string(ans);
        cur=cur+"\n";
        Write(fd, cur.c_str(), cur.size());
    }
    sum /= round;
    LOG_INFO("Raw result has been written to output/exp1-cm.txt");
    LOG_RESULT("Predicted Error Rate: %lf", sum);
    close(fd);

    LOG_INFO("Current Sketch: Count");
    fd=Open("output/exp1-count.txt", O_WRONLY | O_CREAT);
    CountSketch c(1);
    LOG_RESULT("Ground Truth: %lf", c.ground_truth(3, 100000, constraints, stream)[0]);
    sum=0;
    for (int i=0;i<round;i++)
    {
        double ans=c.simulate(3, 100000, constraints, stream, false, false, false)[0];
        sum += ans;
        string cur=to_string(ans);
        cur=cur+"\n";
        Write(fd, cur.c_str(), cur.size());
    }
    sum /= round;
    LOG_INFO("Raw result has been written to output/exp1-count.txt");
    LOG_RESULT("Predicted Error Rate: %lf", sum);
    close(fd);

    LOG_INFO("Current Sketch: Tower");
    fd=Open("output/exp1-tower.txt", O_WRONLY | O_CREAT);
    TowerSketch t(1);
    LOG_RESULT("Ground Truth: %lf", t.ground_truth(3, 100000, constraints, stream)[0]);
    sum=0;
    for (int i=0;i<round;i++)
    {
        double ans=t.simulate(3, 100000, constraints, stream, false, false, false)[0];
        sum += ans;
        string cur=to_string(ans);
        cur=cur+"\n";
        Write(fd, cur.c_str(), cur.size());
    }
    sum /= round;
    LOG_INFO("Raw result has been written to output/exp1-tower.txt");
    LOG_RESULT("Predicted Error Rate: %lf", sum);
    close(fd);

    LOG_INFO("Current Sketch: NitroCM");
    fd=Open("output/exp1-nitrocm.txt", O_WRONLY | O_CREAT);
    NitroCMSketch nitro(1, 0.1);
    LOG_RESULT("Ground Truth: %lf", nitro.ground_truth(3, 200000, constraints, stream)[0]);
    sum=0;
    for (int i=0;i<round;i++)
    {
        double ans=nitro.simulate(3, 200000, constraints, stream, false, false, false)[0];
        sum += ans;
        string cur=to_string(ans);
        cur=cur+"\n";
        Write(fd, cur.c_str(), cur.size());
    }
    sum /= round;
    LOG_INFO("Raw result has been written to output/exp1-nitrocm.txt");
    LOG_RESULT("Predicted Error Rate: %lf", sum);
    close(fd);
}

void experiment_2()
{
    LOG_INFO("Experiment #2: Experiments on Reusing simulation results");
    const char input[] = "input/demand-exp234.txt";
    yyin = fopen(input, "r");
    assert(yyin);
    std::unique_ptr<BaseAst> ast;
    auto ret = yyparse(ast);
    std::cout.flush();
    assert(!ret);
    LOG_INFO("Demands:");
    ast->Print();
    auto demands=reinterpret_cast<ASTDemands*>(ast.get())->Dump();

    TP start, finish;
    StreamGen stream;
    stream.init("dataset/data.dat", 4);

    for (auto& demand : demands)
    {
        LOG_INFO("Current demand: %s", demand.name.c_str());
        LOG_INFO("Without Reuse:");
        start=now();
        search(stream, demand, false, false, false);
        finish=now();
        double t = std::chrono::duration_cast< std::chrono::duration<double, std::ratio<1, 1> > >(finish - start).count();
        LOG_RESULT("time: %lf", t);
        LOG_INFO("With Reuse:");
        start=now();
        search(stream, demand, false, false, true);
        finish=now();
        t = std::chrono::duration_cast< std::chrono::duration<double, std::ratio<1, 1> > >(finish - start).count();
        LOG_RESULT("time: %lf", t);
    }
}

void experiment_3()
{
    LOG_INFO("Experiment #3: Experiments on the Optimization");
    const char input[] = "input/demand-exp234.txt";
    yyin = fopen(input, "r");
    assert(yyin);
    std::unique_ptr<BaseAst> ast;
    auto ret = yyparse(ast);
    std::cout.flush();
    assert(!ret);
    LOG_INFO("Demands:");
    ast->Print();
    auto demands=reinterpret_cast<ASTDemands*>(ast.get())->Dump();

    TP start, finish;
    StreamGen stream;
    stream.init("dataset/data.dat", 4);
    TINY_THRESHOLD=50;
    stream.trunc_stream_init(TINY_THRESHOLD);
    BIG_PERCENT=double(stream.TOTAL_BIG)/stream.TOTAL_FLOWS;

    for (auto& demand : demands)
    {
        LOG_INFO("Current demand: %s", demand.name.c_str());
        LOG_INFO("Without Optimization:");
        start=now();
        search(stream, demand, false, false, true);
        finish=now();
        double t = std::chrono::duration_cast< std::chrono::duration<double, std::ratio<1, 1> > >(finish - start).count();
        LOG_RESULT("time: %lf", t);
        LOG_INFO("With Optimization:");
        start=now();
        search(stream, demand, false, true, true);
        finish=now();
        t = std::chrono::duration_cast< std::chrono::duration<double, std::ratio<1, 1> > >(finish - start).count();
        LOG_RESULT("time: %lf", t);
    }
}

void experiment_4()
{
    LOG_INFO("Experiment #4: Experiments on Configurations");
    const char input[] = "input/demand-exp234.txt";
    yyin = fopen(input, "r");
    assert(yyin);
    std::unique_ptr<BaseAst> ast;
    auto ret = yyparse(ast);
    std::cout.flush();
    assert(!ret);
    LOG_INFO("Demands:");
    ast->Print();
    auto demands=reinterpret_cast<ASTDemands*>(ast.get())->Dump();

    TP start, finish;
    StreamGen stream;
    stream.init("dataset/data.dat", 4);
    TINY_THRESHOLD=50;
    stream.trunc_stream_init(TINY_THRESHOLD);
    BIG_PERCENT=double(stream.TOTAL_BIG)/stream.TOTAL_FLOWS;

    for (auto& demand : demands)
    {
        LOG_INFO("Current demand: %s", demand.name.c_str());
        LOG_INFO("Theory-based Configuration:");
        StrawmanSearch(stream, demand);
        LOG_INFO("Baseline:");
        start=now();
        search(stream, demand, true, false, true);
        finish=now();
        double t = std::chrono::duration_cast< std::chrono::duration<double, std::ratio<1, 1> > >(finish - start).count();
        LOG_RESULT("time: %lf", t);
        LOG_INFO("Ours:");
        start=now();
        search(stream, demand, false, true, true);
        finish=now();
        t = std::chrono::duration_cast< std::chrono::duration<double, std::ratio<1, 1> > >(finish - start).count();
        LOG_RESULT("time: %lf", t);
    }
}

void experiment_5()
{
    LOG_INFO("Experiment #5: Experiments on Sketch Re-configuration");
    LOG_INFO("NOTICE: To run this experiment, you need to place 60 data files named 130000.dat to 135900.dat in the dataset directory!");
    yyin = fopen("input/demand-exp5.txt", "r");
    assert(yyin);
    std::unique_ptr<BaseAst> ast1;
    auto ret1 = yyparse(ast1);
    std::cout.flush();
    assert(!ret1);
    ast1->Print();
    auto demands=reinterpret_cast<ASTDemands*>(ast1.get())->Dump();
    if (demands.size()!=1)
    {
        LOG_ERROR("demands.size=%lu, incorrect!", demands.size());
        exit(-1);
    }
    auto constraints=demands[0].constraints;
    SKETCH::CMSketch* cm[2];
    cm[0]=new SKETCH::CMSketch(4, 85000);
    cm[1]=new SKETCH::CMSketch(4, 85000);

    for (int i=1300; i<=1359; i++)
    {
        TP start, finish;
        TP stt, fih;
        TP ss, ff;
        StreamGen stream;
        stream.init("dataset/"+to_string(i*100)+".dat", 21);
        int use=-1, calc=-1;
        use = i%2;
        if (i>1300)
            calc=1-use;
        
        // insert data
        LOG_INFO("using cm(%d, %d)", cm[use]->nrows_, cm[use]->len_);
        start=now();
        for (int i=0;i<stream.TOTAL_PACKETS;i++)
        {
            cm[use]->insert(stream.raw_data[i]);
        }
        finish=now();
        double t = std::chrono::duration_cast< std::chrono::duration<double, std::ratio<1, 1> > >(finish - start).count();
        LOG_RESULT("thp: %lf", stream.TOTAL_PACKETS/(1e+6 * t));

        vector<double> err(constraints.size(), 0);
        for (auto& it : stream.counter)
        {
            int curst=cm[use]->query(it.first);
            for (int k=0;k<constraints.size();k++)
            {
                if (curst-it.second > constraints[k].err)
                    err[k]++;
            }
        }
        double nflows=stream.getTotalFlows();
        for (int k=0;k<constraints.size();k++)
        {
            LOG_RESULT("err[%d]=%lf", k, err[k]/nflows);
        }

        // analyze previous sketch
        if (calc!=-1)
        {
            start=now();
            EMFSD mrac;
            mrac.set_counters(cm[calc]->len_, cm[calc]->nt_[0]);
            mrac.next_epoch();
            mrac.next_epoch();
            mrac.next_epoch();
            ff=now();
            vector<double> dist=mrac.ns;
            StreamGen prevs;
            prevs.init(dist);
            TINY_THRESHOLD=50;
            prevs.trunc_stream_init(TINY_THRESHOLD);
            BIG_PERCENT=double(prevs.TOTAL_BIG)/prevs.TOTAL_FLOWS;

            stt=now();
            auto conf = search(prevs, demands[0], false, true, true);
            fih=now();

            delete cm[calc];
            cm[calc]= new SKETCH::CMSketch(conf.first, conf.second);
            finish=now();

            double t0 = std::chrono::duration_cast< std::chrono::duration<double, std::ratio<1, 1> > >(finish - start).count();
            LOG_RESULT("t0 = %lf", t0);
            double t1 = std::chrono::duration_cast< std::chrono::duration<double, std::ratio<1, 1> > >(ff - start).count();
            LOG_RESULT("t1 = %lf", t1);
            double t2 = std::chrono::duration_cast< std::chrono::duration<double, std::ratio<1, 1> > >(fih - stt).count();
            LOG_RESULT("t2 = %lf", t2);
        }
        else
            LOG_INFO("No analyze task");
        
    }
}
