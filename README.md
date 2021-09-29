---
R code for the paper entitled Lower circulating lymphocyte count predicts APOE ε4-related cognitive decline in Parkinson’s disease
---
This is the original R code for the manuscript whose DOI is "10.1002/mds.28799". For any question, you can contact Kazuto Tsukita (kazusan@kuhp.kyoto-u.ac.jp).

<br />  

##**Preparation for analysis**  
We downloaded all data to "~/download/" from PPMI website. We used following libraries for subsequent analyses.

    suppressMessages(library(tidyverse)); library(stringr)  %>% suppressMessages(); library(gtsummary)  %>% suppressMessages(); library(lubridate)  %>%           suppressMessages();library(lme4)  %>% suppressMessages();library(plyr)  %>% suppressMessages();library(Amelia) %>% suppressMessages()
    options(gtsummary.pvalue_fun = function(x) style_pvalue(x, digits = 2))
  <br />  
  
##**Identifying PD patients**  
We first identified PD patients from "Primary_Diagnosis.csv".  
In the "Primary_Diagnosis.csv", Idiopathic PD was coded with "01" in the "PRIMDIAG" column.  
  
    PrimaryDiagnosis_PD<-readr::read_csv("~/download/Primary_Diagnosis.csv") %>% dplyr::filter(PRIMDIAG=="01") %>% suppressMessages() %>% suppressWarnings()
    PrimaryDiagnosis_PD<-PrimaryDiagnosis_PD[order(PrimaryDiagnosis_PD$PATNO,PrimaryDiagnosis_PD$EVENT_ID),]

A few patients undergo 2 evaluations with the same "EVENT_ID", and this causes a problem when this dataset were combined with others by "PATNO" and "EVENT_ID"; therefore, we only used the first one.
    c<-c()
    for(i in seq(2,nrow(PrimaryDiagnosis_PD))){
      if((PrimaryDiagnosis_PD$PATNO[i]==PrimaryDiagnosis_PD$PATNO[i-1]) &(PrimaryDiagnosis_PD$EVENT_ID[i]==PrimaryDiagnosis_PD$EVENT_ID[i-1]))(append(c,"NG")->c)
      else(append(c,"OK")->c)
    }
    append("OK",c)->c; check<-c;rm(c,i)

    PrimaryDiagnosis_PD<-cbind(PrimaryDiagnosis_PD,check) %>% dplyr::filter(check=="OK") %>% dplyr::select(-check);rm(check)
  
The number of PD patients enrolled at the baseline in original PPMI study was checked.  
    PrimaryDiagnosis_PD %>% dplyr::filter(EVENT_ID=="SC") %>% with(unique(PATNO)) %>% length()
  <br />
  <br />
##**Identifying PD patients with the information related to "APOE genotype**  
We then identified PD patients with the information related to "APOE genotype"  

    bio<-readr::read_csv("~/download/Current_Biospecimen_Analysis_Results.csv") %>% suppressMessages() 
    ApoE<-bio %>% dplyr::filter(TESTNAME=="ApoE Genotype") %>% dplyr::select(PATNO,TESTVALUE)
    colnames(ApoE)[colnames(ApoE)=="TESTVALUE"]<-"ApoE_Genotype"

    PrimaryDiagnosis_PD_APOE <- PrimaryDiagnosis_PD %>% dplyr::inner_join(ApoE, by="PATNO") %>% dplyr::select(PATNO,EVENT_ID,ApoE_Genotype) ;     rm(PrimaryDiagnosis_PD,bio)

The number of PD patients with the information related to "APOE genotype" was checked.  

    PrimaryDiagnosis_PD_APOE %>% dplyr::filter(EVENT_ID=="SC") %>% with(unique(PATNO)) %>% length()
  <br />
  <br />
  
##**Identifying PD patients both with "APOE genotype" and "baseline MoCA score"**  
We then identified PD patients with the information related to "APOE genotype" and baseline MoCA score.  

    MOCA<-read_csv("~/download/Montreal_Cognitive_Assessment__MoCA_.csv") %>% dplyr::select(PATNO,EVENT_ID,ORIG_ENTRY,MCATOT) %>% suppressMessages() %>%     suppressWarnings()
    colnames(MOCA)<-c("PATNO","EVENT_ID","ORIG_ENTRY","MOCA_Total_score")
    MOCA$ORIG_ENTRY<-as.POSIXct(paste(substr(MOCA$ORIG_ENTRY,4,7),substr(MOCA$ORIG_ENTRY,1,2),01,sep="-"))
    PrimaryDiagnosis_PD_APOE_base<-PrimaryDiagnosis_PD_APOE %>% dplyr::filter(EVENT_ID=="SC")
    MOCA_PD_APOE<-MOCA %>% dplyr::filter(is.element(PATNO,PrimaryDiagnosis_PD_APOE_base$PATNO)) %>% dplyr::inner_join(ApoE, by="PATNO")
    rm(MOCA,PrimaryDiagnosis_PD_APOE,ApoE,PrimaryDiagnosis_PD_APOE_base)

The number of PD patients enrolled with the information related to "APOE genotype"and baseline MoCA score was checked.  

    MOCA_PD_APOE %>% dplyr::filter(EVENT_ID=="SC") %>% with(unique(PATNO)) %>% length()
  <br />
  <br />
  
##**Identifying PD patients with "APOE genotype", "baseline MoCA score", and "baseline Lymphocyte count"**  
We finally identified PD patients with information related to APOE genotype, baseline MoCA score, and baseline Lymphocyte count for this study.
```{r}
blood<-readr::read_csv("~/download/Blood_Chemistry___Hematology.csv") %>% suppressMessages()
Lymphocytes <- blood %>% dplyr::filter(LTSTNAME=="Lymphocytes") %>% dplyr::filter(is.na(LSIRES)==F)
Lymphocyte_BL <- Lymphocytes %>% dplyr::filter(EVENT_ID=="SC") %>% dplyr::select(PATNO,LUSRES)
rm(blood, Lymphocytes)
```
<br />
A few patients undergo 2 measurement with the same "EVENT_ID", and this causes a problem when this dataset were combined with others by "PATNO" and "EVENT_ID"; therefore, we only used the first one.  

    Lymphocyte_BL_1<-c()
    for(i in unique(Lymphocyte_BL$PATNO)){
      Lymphocyte_BL %>% dplyr::filter(PATNO==i)->a1;
      Lymphocyte_BL_1<-bind_rows(Lymphocyte_BL_1,a1[1,])
    }
    Lymphocyte_BL<-Lymphocyte_BL_1
    rm(Lymphocyte_BL_1,i,a1)

    Lymphocyte_BL$LUSRES<-as.numeric(Lymphocyte_BL$LUSRES)
    colnames(Lymphocyte_BL)[colnames(Lymphocyte_BL)=="LUSRES"]<-"BL_Lymphocytes"

    MOCA_PD_APOE_LymphoBL<-MOCA_PD_APOE %>% dplyr::inner_join(Lymphocyte_BL,by="PATNO")
    rm(Lymphocyte_BL,MOCA_PD_APOE)

<br />  
The number of PD patients with "APOE genotype", "baseline MoCA score", and "baseline Lymphocyte count" was checked.  
    MOCA_PD_APOE_LymphoBL %>% dplyr::filter(EVENT_ID=="SC") %>% with(unique(PATNO)) %>% length()  
<br />  
<br />  

##**Incorporating other cognitive tests into the dataset**  
We also used the result of HVLT-R, JLO, and SDMT. These test were first conducted at EVENT_ID of "BL", although MoCA tests were first conducted at EVENT_ID of "SC". These two EVENT_IDs were separated only by 45 days. Therefore, at first, we changed the "SC" to "BL" in the dataset, followed by incorporating HVLT-R, JLO, and SDMT into the dataset  

    MOCA_PD_APOE_LymphoBL$EVENT_ID<-replace(MOCA_PD_APOE_LymphoBL$EVENT_ID,which(MOCA_PD_APOE_LymphoBL$EVENT_ID=="SC"),"BL")

    HVLT<-read_csv("~/download/Hopkins_Verbal_Learning_Test.csv") %>% suppressMessages() %>% suppressWarnings()
    HVLT<-HVLT %>% dplyr::select("PATNO", "EVENT_ID","DVT_DELAYED_RECALL"); colnames(HVLT)<-c("PATNO", "EVENT_ID","HVLT_DELAYED_RECALL_T")
    MOCA_PD_APOE_LymphoBL_HVLT<-MOCA_PD_APOE_LymphoBL %>% dplyr::left_join(HVLT,by = c("PATNO", "EVENT_ID" ))

    JLO<-read_csv("~/download/Benton_Judgment_of_Line_Orientation.csv") %>% suppressMessages()
    JLO<-JLO %>% dplyr::select("PATNO", "EVENT_ID","PAG_NAME","JLO_TOTCALC")
    JLO_1<-JLO %>% dplyr::filter (PAG_NAME=="BENTONOD")
    JLO_2<-JLO %>% dplyr::filter (PAG_NAME=="BENTONEV")
    JLO_3<-dplyr::inner_join(JLO_1,JLO_2,by = c("PATNO", "EVENT_ID"))
    JLO_3<-dplyr::setdiff(JLO_2[,1:2],JLO_3[,1:2],by = c("PATNO", "EVENT_ID"))
    JLO<-dplyr::left_join(JLO_3,JLO_2,by = c("PATNO", "EVENT_ID")) %>% rbind(JLO_1);rm(JLO_1,JLO_2,JLO_3)
    JLO<-JLO %>% dplyr::select("PATNO", "EVENT_ID","JLO_TOTCALC")
    colnames(JLO)<-c("PATNO", "EVENT_ID","JLO_RAW")
    MOCA_PD_APOE_LymphoBL_HVLT_JLO<-MOCA_PD_APOE_LymphoBL_HVLT %>% dplyr::left_join(JLO,by = c("PATNO", "EVENT_ID"))

    SDMT<-read_csv("~/download/Symbol_Digit_Modalities.csv") %>% suppressMessages()
    SDMT<-SDMT %>% dplyr::select("PATNO", "EVENT_ID","SDMTOTAL")
    colnames(SDMT)<-c("PATNO", "EVENT_ID","SDMT_correct")
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT<-MOCA_PD_APOE_LymphoBL_HVLT_JLO %>% dplyr::left_join(SDMT,by = c("PATNO", "EVENT_ID" ))

    UPSIT<-readr::read_csv("~/download/University_of_Pennsylvania_Smell_ID_Test.csv") %>% dplyr::mutate(UPSIT_total=UPSITBK1+UPSITBK2+UPSITBK3+UPSITBK4) %>% dplyr::select(PATNO,EVENT_ID,UPSIT_total) %>% suppressMessages()

    RBDSQ<-read_csv("~/download/REM_Sleep_Disorder_Questionnaire.csv") %>% suppressMessages()
    RBDSQ<-RBDSQ %>% mutate(RBDSQ_last=if_else((STROKE==0)&(HETRA+PARKISM==0)&(RLS==0)&(NARCLPSY==0)&(DEPRS==0)&(EPILEPSY==0)&(BRNINFM==0)&(CNSOTH==0),0,1))
    RBDSQ<-RBDSQ %>% mutate(RBDSQ_total_score=(PTCGBOTH+DRMVIVID+DRMAGRAC+DRMNOCTB+SLPLMBMV+SLPINJUR+DRMVERBL+DRMFIGHT+DRMUMV+DRMOBJFL+MVAWAKEN+DRMREMEM+SLPDSTRB+RBDSQ_last)) %>% dplyr::select(PATNO,EVENT_ID,RBDSQ_total_score)

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT %>% dplyr::left_join(UPSIT,by = c("PATNO", "EVENT_ID" )) %>% dplyr::left_join(RBDSQ,by = c("PATNO", "EVENT_ID" ))

    rm(HVLT,JLO,SDMT,MOCA_PD_APOE_LymphoBL,MOCA_PD_APOE_LymphoBL_HVLT,MOCA_PD_APOE_LymphoBL_HVLT_JLO,RBDSQ,UPSIT)

<br />  

It is very difficult to interepret the results obtained by unschedulled visits. Therefore, we excluded these data from the dataset.  

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$Elapsed_months<-as.numeric(plyr::mapvalues(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$EVENT_ID, from= c("SC","BL","V01", "V02","V03","V04","V05", "V06", "V07","V08","V09", "V10", "V11","V12","V13","V14", "V15","V16"), to = c(0,0,3,6,9,12,18,24,30,36,42,48,54,60,72,84,96,108))) %>% suppressWarnings()  %>% suppressMessages()
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT<- MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT %>% dplyr::filter(is.na(Elapsed_months)==F)

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT[order(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$PATNO,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$EVENT_ID),]

<br />  
<br />  

##**Adding background characteristics to the dataset**  
We fist added Age and sex.  

    screening_Demographic<-read_csv("~/download/Screening___Demographics.csv") %>% filter(is.element(PATNO,unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$PATNO))) %>% suppressMessages()
    screening_Demographic$BIRTHDT<-as.numeric(screening_Demographic$BIRTHDT)
    screening_Demographic<-screening_Demographic%>%mutate(Gender=if_else(GENDER==2,"Male","Female"))
    screening_Demographic<-screening_Demographic[order(screening_Demographic$PATNO,screening_Demographic$EVENT_ID),]

    c<-c()
    for(i in unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$PATNO)){
      screening_Demographic %>% dplyr::filter(PATNO==i) ->a;
      a$BIRTHDT-> a;
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT %>% dplyr::filter(PATNO==i) ->b;
      as.numeric(format(b$ORIG_ENTRY,"%Y"))->b;
      append(c,b-a)->c
    }
    Age<-c

    c<-c()
    for(i in unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$PATNO)){
      screening_Demographic %>% dplyr::filter(PATNO==i) ->a;
      a$Gender-> a;
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT %>% dplyr::filter(PATNO==i) ->b;
      append(c,rep(a,nrow(b)))->c
    }
    Sex<-c

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT<-cbind(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT[,1:2],Age,Sex,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT[,3:ncol(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT)]) %>% as_tibble()
    rm(Age,Sex,a,c,i,b,screening_Demographic)

<br />  
<br />  

Then, we added disease duration.  

    PD_features<-read_csv("~/download/PD_Features.csv") %>% filter(is.element(PATNO,unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$PATNO))) %>% filter(EVENT_ID=="SC") %>% suppressMessages()
    PD_features$PDDXDT<-as.POSIXct(paste(substr(PD_features$PDDXDT,4,7),substr(PD_features$PDDXDT,1,2),01,sep="-"))
    PD_features<-PD_features[order(PD_features$PATNO,PD_features$EVENT_ID),]

    c<-c()
    for(i in unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$PATNO)){
      PD_features %>% dplyr::filter(PATNO==i) ->a;
      as.numeric(format(a$PDDXDT,"%Y"))-> a;
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT %>% dplyr::filter(PATNO==i) ->b;
      as.numeric(format(b$ORIG_ENTRY,"%Y"))->b;
      append(c,b-a)->c
    }
    Disease_duration_diagnosis<-c;rm(a,b,c,i)

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT<-cbind(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT[,1:4],Disease_duration_diagnosis,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT[,5:ncol(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT)]) %>% as_tibble()
    rm(Disease_duration_diagnosis,PD_features)
    
<br />  

We also added LEDD.

    CONMED<-read_csv("~/download/Concomitant_Medications.csv") %>% filter(PD_MOTOR_MED==1) %>% suppressMessages() %>% suppressWarnings()
    CONMED<-CONMED[is.na(CONMED$ORIG_ENTRY)==F,]
    CONMED$ORIG_ENTRY<-as.POSIXct(paste(substr(CONMED$ORIG_ENTRY,4,7),substr(CONMED$ORIG_ENTRY,1,2),01,sep="-"))
    CONMED<-CONMED %>% dplyr::select(PATNO,LEDD,ORIG_ENTRY);CONMED<-CONMED[order(CONMED$PATNO,CONMED$ORIG_ENTRY),]
    CONMED$LEDD<-as.numeric(CONMED$LEDD) %>% suppressWarnings() %>% suppressMessages()

    c<-c()
    for(i in unique(CONMED$PATNO)){
      CONMED %>% filter(PATNO==i)->a;
      for(s in unique(a$ORIG_ENTRY)){
        a %>% filter(ORIG_ENTRY==s)->b;
        sum(b$LEDD,na.rm=T)->b;
        append(c,b)->c
      }
    }
    LEDD_sum<-c

    c<-c()
    for(i in unique(CONMED$PATNO)){
      CONMED %>% filter(PATNO==i)->a;
      for(s in seq(1,nrow(a))){
        if(nrow(a)==s |(a$ORIG_ENTRY[s]!=a$ORIG_ENTRY[s+1]))(append(c,"OK")->c)else(append(c,"NG")->c)
      }
    }
    check<-c;CONMED<-cbind(CONMED,check) %>% filter(check=="OK");rm(c,check,b,i,s)
    CONMED<-cbind(CONMED[colnames(CONMED)=="PATNO"],CONMED[colnames(CONMED)=="ORIG_ENTRY"],LEDD_sum) %>% as_tibble()
    colnames(CONMED)<-c("PATNO","ORIG_ENTRY_LEDD","LEDD_sum")

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT %>% dplyr::left_join(CONMED,by=c("PATNO"))
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_2<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1 %>% filter(is.na(ORIG_ENTRY_LEDD))
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1 %>% filter(is.na(ORIG_ENTRY_LEDD)==F)

    c<-c()
    for(i in unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1$PATNO)){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1 %>% filter(PATNO==i)->a;
      for(s in unique(a$EVENT_ID)){
        a %>% filter(EVENT_ID==s) ->b;
        for(t in seq(1,nrow(b))){
          if((b$ORIG_ENTRY[t]-b$ORIG_ENTRY_LEDD[t])>0)(append(c,"OK")->c)else(append(c,"NG")->c)
        }
      }
    }
    check_1<-c
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1 %>% cbind(check_1) %>% as_tibble()

    c<-c()
    for(i in unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1$PATNO)){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1 %>% filter(PATNO==i)->a;
      for(s in unique(a$EVENT_ID)){
        a %>% filter(EVENT_ID==s) ->b;
        if(sum(b$check_1=="NG")==nrow(b)){b$LEDD_sum<-rep(0,nrow(b));bind_rows(c,b[1,])->c}else(
          bind_rows(c,b[b$check_1=="OK",][nrow(b[b$check_1=="OK",]),])->c
        )
      }
    }
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1<-c
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT<-rbind(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1[,1:(ncol(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1)-1)],MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_2) %>% as_tibble()
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT[order(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$PATNO,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$EVENT_ID),]
    rm(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_1,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_2,check_1,i,LEDD_sum,s,t,a,b,c,CONMED)

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$LEDD_sum<-replace(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$LEDD_sum,is.na(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$LEDD_sum),0)
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT[order(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$PATNO,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$EVENT_ID),]

<br />  
<br />  

Finally, we added UPDRS part3 score.

    UPDRS_part3<-read_csv("~/download/MDS_UPDRS_Part_III.csv") %>% suppressMessages() %>% suppressWarnings()
    UPDRS_part3<-UPDRS_part3 %>% dplyr::filter(is.element(PATNO,unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT$PATNO)))
    UPDRS_part3<-UPDRS_part3[order(UPDRS_part3$PATNO,UPDRS_part3$EVENT_ID),]
    UPDRS_part3<-UPDRS_part3 %>% dplyr::filter((PAG_NAME=="NUPDRS3")|(PAG_NAME=="NUPDR3OF")) %>% dplyr::filter((is.na(ON_OFF_DOSE)) | (ON_OFF_DOSE==1)) %>% dplyr::select(PATNO,EVENT_ID,NP3SPCH,NP3FACXP,NP3RIGN,NP3RIGRU,NP3RIGLU,PN3RIGRL,NP3RIGLL,NP3FTAPR,NP3FTAPL,NP3HMOVR,NP3HMOVL,NP3PRSPR,NP3PRSPL,NP3TTAPR,NP3TTAPL,NP3LGAGR,NP3LGAGL,NP3RISNG,NP3GAIT,NP3FRZGT,NP3PSTBL,NP3POSTR,NP3BRADY,NP3PTRMR,NP3PTRML,NP3KTRMR,NP3KTRML,NP3RTARU,NP3RTALU,NP3RTARL,NP3RTALL,NP3RTALJ,NP3RTCON,NHY)
    UPDRS_part3<-UPDRS_part3 %>% mutate(sum_UPDRS_part3=NP3SPCH+NP3FACXP+NP3RIGN+NP3RIGRU+NP3RIGLU+PN3RIGRL+NP3RIGLL+NP3FTAPR+NP3FTAPL+NP3HMOVR+NP3HMOVL+NP3PRSPR+NP3PRSPL+NP3TTAPR+NP3TTAPL+NP3LGAGR+NP3LGAGL+NP3RISNG+NP3GAIT+NP3FRZGT+NP3PSTBL+NP3POSTR+NP3BRADY+NP3PTRMR+NP3PTRML+NP3KTRMR+NP3KTRML+NP3RTARU+NP3RTALU+NP3RTARL+NP3RTALL+NP3RTALJ+NP3RTCON) %>% dplyr::select(PATNO,EVENT_ID,sum_UPDRS_part3,NHY)

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT %>% dplyr::left_join(UPDRS_part3,by = c("PATNO", "EVENT_ID" )) 
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS[order(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS$PATNO,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS$EVENT_ID),]

    rm(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT,UPDRS_part3)

<br />  
<br />  

##**Follow-up years**  
We determined follow-up years at first. The table entitled "percentage" corresponds to Supplementary Table 1.  

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS$Elapsed_years <- MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS$Elapsed_months/12

    follow.up<-c()
    for(i in unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS$PATNO)){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS %>% filter(PATNO==i) %>% dplyr::select(PATNO,Elapsed_years)->a1;
      follow.up<-bind_rows(follow.up,a1[nrow(a1),])
    }

    percentage<-c()
    for(i in seq(0,max(follow.up$Elapsed_years))){
      sum(follow.up$Elapsed_years>=i)/nrow(follow.up)->a1;
      percentage<-bind_rows(percentage,data.frame(years=i,number=sum(follow.up$Elapsed_years>=i),percentage=a1))
    }


Then, we calculated acceptable follow-up years.  


    OK<-c()
    for(i in seq(1,nrow(percentage))){
      if(percentage$percentage[i]>=0.8)(OK<-dplyr::bind_rows(OK,percentage[i,]))
    }

    max(OK$years)

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_all<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS %>% dplyr::filter(Elapsed_years<=max(OK$years))
    rm(OK,a1,i,follow.up)


<br />  
<br />  

##**Descripting of demographics of included patients**  
First, we devided the dataset based on the APOE status.  

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 <- MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS %>% dplyr::filter(stringr::str_detect(ApoE_Genotype,"e4")) %>% dplyr::mutate(ApoE_e4_allele="APOE ε4 carrier")
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4 <- MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS %>% dplyr::filter(stringr::str_detect(ApoE_Genotype,"e4")==F) %>% dplyr::mutate(ApoE_e4_allele="Non-carrier")
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS <- dplyr::bind_rows(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4)


Then, we created the table. By this code, we get the information on which the Supplementary Table 2 is based.  
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS <- MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS %>% dplyr::select("PATNO","EVENT_ID","ApoE_Genotype","Age","Sex","Disease_duration_diagnosis","NHY","LEDD_sum","BL_Lymphocytes","MOCA_Total_score","HVLT_DELAYED_RECALL_T","JLO_RAW","SDMT_correct","sum_UPDRS_part3","UPSIT_total","RBDSQ_total_score","ORIG_ENTRY","ORIG_ENTRY_LEDD","Elapsed_months","Elapsed_years","ApoE_e4_allele")

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS$ApoE_Genotype<-factor(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS$ApoE_Genotype,levels=c("e2/e2","e3/e2","e2/e4","e3/e3","e4/e3","e4/e4"))

    Table_bl<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS[,(str_detect(colnames(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS),"ORIG_ENTRY")==F)] %>% dplyr::filter(Elapsed_years==0) %>% dplyr::select(-PATNO,-EVENT_ID,-Elapsed_months,-Elapsed_years) %>%
    gtsummary::tbl_summary(by = `ApoE_e4_allele`,digits = list(all_continuous() ~ 1,all_categorical() ~ 1),type=c(Disease_duration_diagnosis,JLO_RAW,NHY,LEDD_sum) ~"continuous",value=Sex~"Male",label=list(Age~"Age (years)",Sex~"Sex, Male (%)",Disease_duration_diagnosis~"Disease duration from the diagnosis (years)",NHY~"Hoen-Yahr stage",BL_Lymphocytes~"Baseline lymphocyte count (× 1000/μL)",MOCA_Total_score~"MoCA Total score",HVLT_DELAYED_RECALL_T~"HVLT-R delayed recall T score",JLO_RAW~"JLO total score",SDMT_correct~"SDMT total score",sum_UPDRS_part3~"MDS-UPDRS part3 off score",LEDD_sum~"Levodopa equivalent daily dose (mg)",UPSIT_total~"UPSIT total score",RBDSQ_total_score~"RBDSQ total score")) %>% modify_header(label ~ "****") %>% modify_spanning_header(c("stat_1", "stat_2") ~ "**At the baseline**")

    Table_one<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS[,(str_detect(colnames(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS),"ORIG_ENTRY")==F)] %>% dplyr::filter(Elapsed_years==1) %>% dplyr::select(-PATNO,-EVENT_ID,-Elapsed_months,-Elapsed_years) %>% gtsummary::tbl_summary(by = `ApoE_e4_allele`,digits = list(all_continuous() ~ 1,all_categorical() ~ 1),type=c(Disease_duration_diagnosis,JLO_RAW,NHY) ~"continuous",value=Sex~"Male",label=list(Age~"Age (years)",Sex~"Sex, Male (%)",Disease_duration_diagnosis~"Disease duration from the diagnosis (years)",NHY~"Hoen-Yahr stage",BL_Lymphocytes~"Baseline lymphocyte count (× 1000/μL)",MOCA_Total_score~"MoCA Total score",HVLT_DELAYED_RECALL_T~"HVLT-R delayed recall T score",JLO_RAW~"JLO total score",SDMT_correct~"SDMT total score",sum_UPDRS_part3~"MDS-UPDRS part3 off score",LEDD_sum~"Levodopa equivalent daily dose (mg)",UPSIT_total~"UPSIT total score",RBDSQ_total_score~"RBDSQ total score")) %>% modify_header(label ~ "****") %>% modify_spanning_header(c("stat_1", "stat_2") ~ "**1 year follow up**")

    Table_two<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS[,(str_detect(colnames(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS),"ORIG_ENTRY")==F)] %>% dplyr::filter(Elapsed_years==2) %>% dplyr::select(-PATNO,-EVENT_ID,-Elapsed_months,-Elapsed_years) %>% gtsummary::tbl_summary(by = `ApoE_e4_allele`,digits = list(all_continuous() ~ 1,all_categorical() ~ 1),type=c(Disease_duration_diagnosis,JLO_RAW,NHY) ~"continuous",value=Sex~"Male",label=list(Age~"Age (years)",Sex~"Sex, Male (%)",Disease_duration_diagnosis~"Disease duration from the diagnosis (years)",NHY~"Hoen-Yahr stage",BL_Lymphocytes~"Baseline lymphocyte count (× 1000/μL)",MOCA_Total_score~"MoCA Total score",HVLT_DELAYED_RECALL_T~"HVLT-R delayed recall T score",JLO_RAW~"JLO total score",SDMT_correct~"SDMT total score",sum_UPDRS_part3~"MDS-UPDRS part3 off score",LEDD_sum~"Levodopa equivalent daily dose (mg)",UPSIT_total~"UPSIT total score",RBDSQ_total_score~"RBDSQ total score")) %>% modify_header(label ~ "****") %>% modify_spanning_header(c("stat_1", "stat_2") ~ "**2 year follow up**")

    tbl_merge(list(Table_bl,Table_one,Table_two), tab_spanner = c("**At the baseline**", "**1 year follow up**","**2 year follow up**"))
    rm(Table_bl,Table_one,Table_two)

<br />  
<br />  

##**Main Analyses**  
We assessed the interaction effect of lymphocyte count at the baseline on the longitudinal change of global cognitive function. The table entitiled "interaction_effect_table_BL_Lymphocytes_main" corresponds to the main result. As the covariates, we selected age, levodopa equivalent daily dose, disease duration, sex, and baseline scores of UPSIT and RBDSQ.  

    all<-c()
    for (i in unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS$PATNO)){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS %>% filter(PATNO==i) ->a1;
      a1$base_UPSIT_total<-rep(a1$UPSIT_total[1],nrow(a1));
      a1$base_RBDSQ_total_score<-rep(a1$RBDSQ_total_score[1],nrow(a1));
      all<-bind_rows(all,a1);
    }
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS <- all
    rm(all)

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 <- MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS %>% dplyr::filter(stringr::str_detect(ApoE_Genotype,"e4"))
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4 <- MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS %>% dplyr::filter(stringr::str_detect(ApoE_Genotype,"e4")==F)


    cognitive_tests<-c("MOCA_Total_score")

    Partner<-c()
    c<-c()
    for (i in cognitive_tests){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS %>% dplyr::select(i,PATNO,Elapsed_years,BL_Lymphocytes,Age,LEDD_sum,Disease_duration_diagnosis,Sex,base_UPSIT_total,base_RBDSQ_total_score) %>% suppressMessages()->a1;
      colnames(a1)[1]<-"V1";
      a1$V1<-scale(a1$V1,center=T,scale=T);
      a1$BL_Lymphocytes<-scale(a1$BL_Lymphocytes,center=T,scale=T);
      a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
      a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
      a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
      a1$base_UPSIT_total<-scale(a1$base_UPSIT_total,center=T,scale=T);
      a1$base_RBDSQ_total_score<-scale(a1$base_RBDSQ_total_score,center=T,scale=T);
      lme4::lmer(V1~Elapsed_years*BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b
      lme4::lmer(V1~Elapsed_years+BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b.null
      summary(b)->b.summary;confint(b) %>% suppressMessages() ->b.confint;
      b.summary$coefficients[10,]->e;b.confint[12,]->f;b.summary$AICtab[1:2]->A.BIC;
      anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
      c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC)->d;
      bind_rows(c,d) %>% suppressMessages()->c;
      append(Partner,i)->Partner
    }
    colnames(c)[6]<-"p.value"
    cbind(Partner,c) %>% dplyr::mutate(Label="all") ->interaction_effect_table_BL_Lymphocytes

    Partner<-c()
    c<-c()
    for (i in cognitive_tests){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 %>% dplyr::select(i,PATNO,Elapsed_years,BL_Lymphocytes,Age,LEDD_sum,Disease_duration_diagnosis,Sex,base_UPSIT_total,base_RBDSQ_total_score)->a1;
      colnames(a1)[1]<-"V1";
      a1$V1<-scale(a1$V1,center=T,scale=T);
      a1$BL_Lymphocytes<-scale(a1$BL_Lymphocytes,center=T,scale=T);
      a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
      a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
      a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
      a1$base_UPSIT_total<-scale(a1$base_UPSIT_total,center=T,scale=T);
      a1$base_RBDSQ_total_score<-scale(a1$base_RBDSQ_total_score,center=T,scale=T);
      lme4::lmer(V1~Elapsed_years*BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b
      lme4::lmer(V1~Elapsed_years+BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b.null
      summary(b)->b.summary;confint(b) %>% suppressMessages() ->b.confint;
      b.summary$coefficients[10,]->e;b.confint[12,]->f;b.summary$AICtab[1:2]->A.BIC;
      anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
      c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC)->d;
      bind_rows(c,d) %>% suppressMessages()->c;
      append(Partner,i)->Partner
    }
    colnames(c)[6]<-"p.value"
    cbind(Partner,c) %>% as_tibble() %>% dplyr::mutate(Label="e4") ->interaction_effect_table_BL_Lymphocytes_e4 

    Partner<-c()
    c<-c()
    for (i in cognitive_tests){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4 %>% dplyr::select(i,PATNO,Elapsed_years,BL_Lymphocytes,Age,LEDD_sum,Disease_duration_diagnosis,Sex,base_UPSIT_total,base_RBDSQ_total_score)->a1;
      colnames(a1)[1]<-"V1";
      a1$V1<-scale(a1$V1,center=T,scale=T);
      a1$BL_Lymphocytes<-scale(a1$BL_Lymphocytes,center=T,scale=T);
      a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
      a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
      a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
      a1$base_UPSIT_total<-scale(a1$base_UPSIT_total,center=T,scale=T);
      a1$base_RBDSQ_total_score<-scale(a1$base_RBDSQ_total_score,center=T,scale=T);
      lme4::lmer(V1~Elapsed_years*BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b
      lme4::lmer(V1~Elapsed_years+BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b.null
      summary(b)->b.summary;confint(b) %>% suppressMessages() ->b.confint;
      b.summary$coefficients[10,]->e;b.confint[12,]->f;b.summary$AICtab[1:2]->A.BIC;
      anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
      c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC) %>% suppressMessages()->d;
      bind_rows(c,d) %>% suppressMessages()->c;
      append(Partner,i)->Partner
    }
    colnames(c)[6]<-"p.value"
    cbind(Partner,c) %>% as_tibble() %>% dplyr::mutate(Label="note4") ->interaction_effect_table_BL_Lymphocytes_note4

    interaction_effect_table_BL_Lymphocytes<-bind_rows(interaction_effect_table_BL_Lymphocytes,interaction_effect_table_BL_Lymphocytes_e4,interaction_effect_table_BL_Lymphocytes_note4)
    interaction_effect_table_BL_Lymphocytes_main<-interaction_effect_table_BL_Lymphocytes %>% dplyr::mutate(significance=if_else(p.value<=0.05,"Sig","Not"))

    rm(a1,b,b.confint,b.null,b.summary,c,d.summary,interaction_effect_table_BL_Lymphocytes,interaction_effect_table_BL_Lymphocytes_e4,interaction_effect_table_BL_Lymphocytes_note4,A.BIC,cognitive_tests,d,e,f,i,Partner)

<br />  

For visualization, both PD patients with and without APOE ε4 allele were dichotomized using the median of lymphocyte count(please refer the vectors named "median_ApoE4" and "median_notApoE4") at the baseline. After this code, we will get the information on which Table 1 is based.  

    median_ApoE4<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 %>% dplyr::filter(EVENT_ID=="BL") %>% with(median(BL_Lymphocytes))
    median_notApoE4<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4 %>% dplyr::filter(EVENT_ID=="BL") %>% with(median(BL_Lymphocytes))

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4$ApoE_Genotype<-factor(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4$ApoE_Genotype,levels=c("e2/e2","e3/e2","e2/e4","e3/e3","e4/e3","e4/e4"))
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4$ApoE_Genotype<-factor(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4$ApoE_Genotype,levels=c("e2/e2","e3/e2","e2/e4","e3/e3","e4/e3","e4/e4"))

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 %>% dplyr::mutate(LymphoBL_cat=if_else(BL_Lymphocytes>=median_ApoE4,"Higher","Lower")) ->a
    a$LymphoBL_cat <- factor(a$LymphoBL_cat,levels=c("Lower","Higher"))
    Table_1_1 <- a %>% dplyr::filter(Elapsed_years==0) %>% dplyr::select("BL_Lymphocytes","ApoE_Genotype","Age","Sex","Disease_duration_diagnosis","NHY","LEDD_sum","MOCA_Total_score","HVLT_DELAYED_RECALL_T","JLO_RAW","SDMT_correct","sum_UPDRS_part3","LymphoBL_cat","UPSIT_total","RBDSQ_total_score") %>% gtsummary::tbl_summary(by = `LymphoBL_cat`,digits = list(all_continuous() ~ 1,all_categorical()~1),type=list(Disease_duration_diagnosis~"continuous",NHY~"continuous",JLO_RAW~"continuous",LEDD_sum~"continuous"),value=Sex~"Male",label=list(Age~"Age (years)",Sex~"Sex, Male (%)",Disease_duration_diagnosis~"Disease duration from the diagnosis (years)",NHY~"Hoen-Yahr stage",BL_Lymphocytes~"Baseline lymphocyte count (× 1000/μL)",MOCA_Total_score~"MoCA Total score",HVLT_DELAYED_RECALL_T~"HVLT-R delayed recall T score",JLO_RAW~"JLO total score",SDMT_correct~"SDMT total score",sum_UPDRS_part3~"MDS-UPDRS part3 off score",LEDD_sum~"Levodopa equivalent daily dose (mg)",UPSIT_total~"UPSIT total score",RBDSQ_total_score~"RBDSQ total score")) %>% add_p() %>% suppressMessages()

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4 %>% dplyr::mutate(LymphoBL_cat=if_else(BL_Lymphocytes>=median_ApoE4,"Higher","Lower"))->a
    a$LymphoBL_cat <- factor(a$LymphoBL_cat,levels=c("Lower","Higher"))

    Table_1_2 <- a %>% dplyr::filter(Elapsed_years==0) %>% dplyr::select("BL_Lymphocytes","ApoE_Genotype","Age","Sex","Disease_duration_diagnosis","NHY","LEDD_sum","MOCA_Total_score","HVLT_DELAYED_RECALL_T","JLO_RAW","SDMT_correct","sum_UPDRS_part3","LymphoBL_cat","UPSIT_total","RBDSQ_total_score") %>% gtsummary::tbl_summary(by = `LymphoBL_cat`,digits = list(all_continuous() ~ 1,all_categorical()~1),type=list(Disease_duration_diagnosis~"continuous",NHY~"continuous",JLO_RAW~"continuous",c(LEDD_sum,MOCA_Total_score)~"continuous"),value=Sex~"Male",label=list(Age~"Age (years)",Sex~"Sex, Male (%)",Disease_duration_diagnosis~"Disease duration from the diagnosis (years)",NHY~"Hoen-Yahr stage",BL_Lymphocytes~"Baseline lymphocyte count (× 1000/μL)",MOCA_Total_score~"MoCA Total score",HVLT_DELAYED_RECALL_T~"HVLT-R delayed recall T score",JLO_RAW~"JLO total score",SDMT_correct~"SDMT total score",sum_UPDRS_part3~"MDS-UPDRS part3 off score",LEDD_sum~"Levodopa equivalent daily dose (mg)",UPSIT_total~"UPSIT total score",RBDSQ_total_score~"RBDSQ total score")) %>% add_p() %>% suppressMessages()

    tbl_merge(list(Table_1_1,Table_1_2), tab_spanner = c("**APOE ε4 carrier**", "**Non-carrier**")) %>% modify_header(label ~ "****") %>% suppressMessages()
    rm(Table_1_1,Table_1_2)

<br />  

Then, temporal change of MoCA total score was plotted. After this code, you will get the information on which Fig.1A and 1B were created.  

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 %>% dplyr::mutate(LymphoBL_cat=if_else(BL_Lymphocytes>=median_ApoE4,"Upper","Lower")) %>% ggplot(aes(x=Elapsed_years,y=MOCA_Total_score,colour=LymphoBL_cat)) + theme_classic()  + stat_smooth(se=T,fullrange = T,level = 0.95,method=lm,alpha=0.2,size=2)+scale_color_manual(name="Baseline lymphocyte count",values = c("#00BFC4","#F8766D"))+scale_y_continuous(breaks = seq(20,30,1))+coord_cartesian(ylim = c(22,29))+scale_x_continuous(breaks = seq(0,2,1)) + ylab("MoCA Total score") + xlab("Follow-up years") + ggtitle("APOE ε4 allele carriers") + theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),legend.background = element_rect(color = "black", linetype = "dashed"),legend.position = c(0.2, 0.15),legend.title = element_text(size = 10),legend.text = element_text(size = 10))

    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4 %>% dplyr::mutate(LymphoBL_cat=if_else(BL_Lymphocytes>=median_ApoE4,"Upper","Lower")) %>% ggplot(aes(x=Elapsed_years,y=MOCA_Total_score,colour=LymphoBL_cat)) + theme_classic()  + stat_smooth(se=T,fullrange = T,level = 0.95,method=lm,alpha=0.2,size=2)+scale_color_manual(name="Baseline lymphocyte count",values = c("#00BFC4","#F8766D"))+scale_y_continuous(breaks = seq(20,30,1))+coord_cartesian(ylim = c(22,29))+scale_x_continuous(breaks = seq(0,2,1)) + ylab("MoCA Total score") + xlab("Follow-up years") + ggtitle("APOE ε4 allele non-carriers") + theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),legend.background = element_rect(color = "black", linetype = "dashed"),legend.position = c(0.2, 0.15),legend.title = element_text(size = 10),legend.text = element_text(size = 10))

<br />  

We also cheched the sifnificance. After this code, we will get the "interaction_effect_table_BL_Lymphocytes_binary_main" containing the result of the multivariate linear mixed-effects models.  

    Partner<-c()
    c<-c()
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 %>% dplyr::mutate(LymphoBL_cat=if_else(BL_Lymphocytes>=median_ApoE4,"Upper","Lower")) %>% dplyr::select(MOCA_Total_score,PATNO,Elapsed_years,LymphoBL_cat,Age,LEDD_sum,Disease_duration_diagnosis,Sex,base_UPSIT_total,base_RBDSQ_total_score)->a1;
    a1$MOCA_Total_score<-scale(a1$MOCA_Total_score,center=T,scale=T);
    a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
    a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
    a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
    a1$base_UPSIT_total<-scale(a1$base_UPSIT_total,center=T,scale=T);
    a1$base_RBDSQ_total_score<-scale(a1$base_RBDSQ_total_score,center=T,scale=T);
    lme4::lmer(MOCA_Total_score~Elapsed_years*LymphoBL_cat+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b
    lme4::lmer(MOCA_Total_score~Elapsed_years+LymphoBL_cat+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b.null
    summary(b)->b.summary;confint(b) %>% suppressMessages() ->b.confint;
    b.summary$coefficients[10,]->e;b.confint[12,]->f;b.summary$AICtab[1:2]->A.BIC;
    anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
    c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC)->d;
    bind_rows(c,d) %>% suppressMessages() ->c;
    append(Partner,"Baseline Lymphocyte caount Binary")->Partner
    colnames(c)[6]<-"p.value"
    cbind(Partner,c) %>% as_tibble() %>% dplyr::mutate(Label="e4")->interaction_effect_table_BL_Lymphocytes_binary_main_1


    Partner<-c()
    c<-c()
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4 %>% dplyr::mutate(LymphoBL_cat=if_else(BL_Lymphocytes>=median_ApoE4,"Upper","Lower")) %>% dplyr::select(MOCA_Total_score,PATNO,Elapsed_years,LymphoBL_cat,Age,LEDD_sum,Disease_duration_diagnosis,Sex,base_UPSIT_total,base_RBDSQ_total_score)->a1;
    a1$MOCA_Total_score<-scale(a1$MOCA_Total_score,center=T,scale=T);
    a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
    a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
    a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
    a1$base_UPSIT_total<-scale(a1$base_UPSIT_total,center=T,scale=T);
    a1$base_RBDSQ_total_score<-scale(a1$base_RBDSQ_total_score,center=T,scale=T);
    lme4::lmer(MOCA_Total_score~Elapsed_years*LymphoBL_cat+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b
    lme4::lmer(MOCA_Total_score~Elapsed_years+LymphoBL_cat+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b.null
    summary(b)->b.summary;confint(b) %>% suppressMessages() ->b.confint;
    b.summary$coefficients[10,]->e;b.confint[12,]->f;b.summary$AICtab[1:2]->A.BIC;
    anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
    c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC)->d;
    bind_rows(c,d) %>% suppressMessages() ->c;
    append(Partner,"Baseline Lymphocyte caount Binary")->Partner
    colnames(c)[6]<-"p.value"
    cbind(Partner,c) %>% as_tibble() %>% dplyr::mutate(Label="note4")->interaction_effect_table_BL_Lymphocytes_binary_main_2

    interaction_effect_table_BL_Lymphocytes_binary_main <- dplyr::bind_rows(interaction_effect_table_BL_Lymphocytes_binary_main_1,interaction_effect_table_BL_Lymphocytes_binary_main_2)

    rm(interaction_effect_table_BL_Lymphocytes_binary_main_1,interaction_effect_table_BL_Lymphocytes_binary_main_2,a1,b,b.confint,b.null,b.summary,c,d.summary,A.BIC,d,e,f,Partner)

<br />  
<br />  

##**Sub-Analyses**  
We further assessed the interaction effect of lymphocyte count at the baseline on the longitudinal change of subdomains of cognitive function. After doing the code, we will get the Fig.1C. The table entitiled "interaction_effect_table_BL_Lymphocytes_main" corresponds to the result.  

    cognitive_tests<-c("HVLT_DELAYED_RECALL_T","JLO_RAW","SDMT_correct")


    Partner<-c()
    c<-c()
    for (i in cognitive_tests){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 %>% dplyr::select(i,PATNO,Elapsed_years,BL_Lymphocytes,Age,LEDD_sum,Disease_duration_diagnosis,Sex,base_UPSIT_total,base_RBDSQ_total_score)->a1;
      colnames(a1)[1]<-"V1";
      a1$V1<-scale(a1$V1,center=T,scale=T);
      a1$BL_Lymphocytes<-scale(a1$BL_Lymphocytes,center=T,scale=T);
      a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
      a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
      a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
      a1$base_UPSIT_total<-scale(a1$base_UPSIT_total,center=T,scale=T);
      a1$base_RBDSQ_total_score<-scale(a1$base_RBDSQ_total_score,center=T,scale=T);
      lme4::lmer(V1~Elapsed_years*BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b
      lme4::lmer(V1~Elapsed_years+BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b.null
      summary(b)->b.summary;confint(b) %>% suppressMessages()->b.confint;
      b.summary$coefficients[10,]->e;b.confint[12,]->f;b.summary$AICtab[1:2]->A.BIC;
      anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
      c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC)->d;
      bind_rows(c,d) %>% suppressMessages()->c;
      append(Partner,i)->Partner
    }
    colnames(c)[6]<-"p.value"
    cbind(Partner,c) %>% as_tibble() %>% dplyr::mutate(Label="e4") ->interaction_effect_table_BL_Lymphocytes_e4 

    Partner<-c()
    c<-c()
    for (i in cognitive_tests){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4 %>% dplyr::select(i,PATNO,Elapsed_years,BL_Lymphocytes,Age,LEDD_sum,Disease_duration_diagnosis,Sex,base_UPSIT_total,base_RBDSQ_total_score)->a1;
      colnames(a1)[1]<-"V1";
      a1$V1<-scale(a1$V1,center=T,scale=T);
      a1$BL_Lymphocytes<-scale(a1$BL_Lymphocytes,center=T,scale=T);
      a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
      a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
      a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
      a1$base_UPSIT_total<-scale(a1$base_UPSIT_total,center=T,scale=T);
      a1$base_RBDSQ_total_score<-scale(a1$base_RBDSQ_total_score,center=T,scale=T);
      lme4::lmer(V1~Elapsed_years*BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b
      lme4::lmer(V1~Elapsed_years+BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b.null
      summary(b)->b.summary;confint(b) %>% suppressMessages()->b.confint;
      b.summary$coefficients[10,]->e;b.confint[12,]->f;b.summary$AICtab[1:2]->A.BIC;
      anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
      c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC)->d;
      bind_rows(c,d) %>% suppressMessages()->c;
      append(Partner,i)->Partner
    }
    colnames(c)[6]<-"p.value"
    cbind(Partner,c) %>% as_tibble() %>% dplyr::mutate(Label="note4") ->interaction_effect_table_BL_Lymphocytes_note4

    interaction_effect_table_BL_Lymphocytes<-bind_rows(interaction_effect_table_BL_Lymphocytes_e4,interaction_effect_table_BL_Lymphocytes_note4)
    interaction_effect_table_BL_Lymphocytes_subdomains<-interaction_effect_table_BL_Lymphocytes %>% dplyr::mutate(significance=if_else(p.value<=0.05,"Sig","Not"))
    interaction_effect_table_BL_Lymphocytes$Partner<-factor(interaction_effect_table_BL_Lymphocytes$Partner,levels=c("JLO_RAW","SDMT_correct","HVLT_DELAYED_RECALL_T","MOCA_Total_score"))
    interaction_effect_table_BL_Lymphocytes_subdomains %>% dplyr::bind_rows(interaction_effect_table_BL_Lymphocytes_main[2:3,])->a
    a$Partner<-factor(a$Partner,levels=c("JLO_RAW","SDMT_correct","HVLT_DELAYED_RECALL_T","MOCA_Total_score"))
    a$Partner<-plyr::mapvalues(a$Partner,from=c("JLO_RAW","SDMT_correct","HVLT_DELAYED_RECALL_T","MOCA_Total_score"),to=c("JLO total score","SDMT total score","HVLT-R delayed recall T score","MoCA Total score"))
    a$Label<-factor(a$Label,levels=c("e4","note4"))
    a$Label<-plyr::mapvalues(a$Label,from=c("e4","note4"),to=c("APOE ε4 carriers","APOE ε4 non-carriers"))
    a%>% ggplot(aes(x=Label,y=Partner,fill=abs(`t value`))) + theme_classic()+ geom_tile()+scale_fill_gradient(low="white",high="red",name="t value",limits=c(0,3.5))+geom_point(aes(shape=significance))+scale_shape_manual(values = c(32,8),labels=c("","significant",""),name="")+ylab("")+xlab("")

    rm(a, a1,b,b.confint,b.null,b.summary,c,d.summary,interaction_effect_table_BL_Lymphocytes,interaction_effect_table_BL_Lymphocytes_e4,interaction_effect_table_BL_Lymphocytes_note4,A.BIC,cognitive_tests,d,e,f,i,Partner)

<br />  
<br />  

##**Sensitivity analyses(1)**  

We conducted sensitivity analyses using different follow-up periods. After this code, we will get the table entitled "interaction_effect_table_BL_Lymphocytes_sensitivity_analyses" which corresponds to Supplementary Table 4.  


    all<-c()
    for (i in unique(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_all$PATNO)){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_all %>% filter(PATNO==i) ->a1;
      a1$base_UPSIT_total<-rep(a1$UPSIT_total[1],nrow(a1));
      a1$base_RBDSQ_total_score<-rep(a1$RBDSQ_total_score[1],nrow(a1));
      all<-bind_rows(all,a1);
    }
    MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_all <- all
    rm(all)

    interaction_effect_table_BL_Lymphocytes_sensitivity_analyses<-c()
    for (s in seq(3,9)){
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_all %>% dplyr::filter(is.element(PATNO,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4$PATNO)) %>% dplyr::filter(Elapsed_years<=s) %>% dplyr::mutate(ApoE_e4_allele="Yes")
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_all %>% dplyr::filter(is.element(PATNO,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4$PATNO)) %>% dplyr::filter(Elapsed_years<=s) %>% dplyr::mutate(ApoE_e4_allele="No")
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS<-dplyr::bind_rows(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4)
      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS$BL_Lymphocytes<-log(MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS$BL_Lymphocytes)
  
      Partner<-c()
      c<-c()
      for (i in "MOCA_Total_score"){
        MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS %>% dplyr::select(i,PATNO,Elapsed_years,BL_Lymphocytes,Age,LEDD_sum,Disease_duration_diagnosis,Sex,base_UPSIT_total,base_RBDSQ_total_score)->a1;
        colnames(a1)[1]<-"V1";
        a1$V1<-scale(a1$V1,center=T,scale=T);
        a1$BL_Lymphocytes<-scale(a1$BL_Lymphocytes,center=T,scale=T);
        a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
        a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
        a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
        a1$base_UPSIT_total<-scale(a1$base_UPSIT_total,center=T,scale=T);
        a1$base_RBDSQ_total_score<-scale(a1$base_RBDSQ_total_score,center=T,scale=T)
        lme4::lmer(V1~Elapsed_years*BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b
        lme4::lmer(V1~Elapsed_years+BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b.null
        summary(b)->b.summary;confint(b) %>% suppressMessages()->b.confint;
        b.summary$coefficients[10,]->e;b.confint[12,]->f;b.summary$AICtab[1:2]->A.BIC;
        anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
        c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC)->d;
        bind_rows(c,d) %>% suppressMessages()->c;
        append(Partner,i)->Partner
      }
      colnames(c)[6]<-"p.value"
      cbind(Partner,c) %>% dplyr::mutate(Label="all") ->interaction_effect_table_BL_Lymphocytes
  
      Partner<-c()
      c<-c()
      for (i in "MOCA_Total_score"){
        MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 %>% dplyr::select(i,PATNO,Elapsed_years,BL_Lymphocytes,Age,LEDD_sum,Disease_duration_diagnosis,Sex)->a1;
        colnames(a1)[1]<-"V1";
        a1$V1<-scale(a1$V1,center=T,scale=T);
        a1$BL_Lymphocytes<-scale(a1$BL_Lymphocytes,center=T,scale=T);
        a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
        a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
        a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
        lme4::lmer(V1~Elapsed_years*BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+(1|PATNO),a1,REML=FALSE)->b
        lme4::lmer(V1~Elapsed_years+BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+(1|PATNO),a1,REML=FALSE)->b.null
        summary(b)->b.summary;confint(b) %>% suppressMessages()->b.confint;
        b.summary$coefficients[8,]->e;b.confint[10,]->f;b.summary$AICtab[1:2]->A.BIC;
        anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
        c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC)->d;
        bind_rows(c,d) %>% suppressMessages()->c;
        append(Partner,i)->Partner
      }
      colnames(c)[6]<-"p.value"
      cbind(Partner,c) %>% as_tibble() %>% dplyr::mutate(Label="e4") ->interaction_effect_table_BL_Lymphocytes_e4 
  
      Partner<-c()
      c<-c()
      for (i in "MOCA_Total_score"){
        MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4 %>% dplyr::select(i,PATNO,Elapsed_years,BL_Lymphocytes,Age,LEDD_sum,Disease_duration_diagnosis,Sex)->a1;
        colnames(a1)[1]<-"V1";
        a1$V1<-scale(a1$V1,center=T,scale=T);
        a1$BL_Lymphocytes<-scale(a1$BL_Lymphocytes,center=T,scale=T);
        a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
        a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
        a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
        lme4::lmer(V1~Elapsed_years*BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+(1|PATNO),a1,REML=FALSE)->b
        lme4::lmer(V1~Elapsed_years+BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+(1|PATNO),a1,REML=FALSE)->b.null
        summary(b)->b.summary;confint(b) %>% suppressMessages()->b.confint;
        b.summary$coefficients[8,]->e;b.confint[10,]->f;b.summary$AICtab[1:2]->A.BIC;
        anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
        c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC)->d;
        bind_rows(c,d) %>% suppressMessages()->c;
        append(Partner,i)->Partner
      }
      colnames(c)[6]<-"p.value"
      cbind(Partner,c) %>% as_tibble() %>% dplyr::mutate(Label="note4") ->interaction_effect_table_BL_Lymphocytes_note4
  
      interaction_effect_table_BL_Lymphocytes<-bind_rows(interaction_effect_table_BL_Lymphocytes,interaction_effect_table_BL_Lymphocytes_e4,interaction_effect_table_BL_Lymphocytes_note4)
      interaction_effect_table_BL_Lymphocytes<-interaction_effect_table_BL_Lymphocytes %>% dplyr::mutate(significance=if_else(p.value<=0.05,"Sig","Not")) %>% dplyr::mutate(follow_up_years=s)
      interaction_effect_table_BL_Lymphocytes_sensitivity_analyses<-dplyr::bind_rows(interaction_effect_table_BL_Lymphocytes_sensitivity_analyses,interaction_effect_table_BL_Lymphocytes)
    }

    rm(a1,b,b.confint,b.null,b.summary,c,d.summary,interaction_effect_table_BL_Lymphocytes,interaction_effect_table_BL_Lymphocytes_e4,interaction_effect_table_BL_Lymphocytes_note4,A.BIC,d,e,f,i,Partner,s,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS,MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_note4)


<br />  
<br />  

##**Sensitivity analyses(2)**  
Finally, we also conducted sensitivity analyses using multiple imputation methodology. After this code, we will get the table entitled "interaction_effect_table_BL_Lymphocytes_sensitivity_analyses" which corresponds to Supplementary Table 5.  

      MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 <- MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 %>% dplyr::filter(Elapsed_years<=2)
      imp<-MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4 <- MOCA_PD_APOE_LymphoBL_HVLT_JLO_SDMT_UPDRS_e4  %>% dplyr::select(MOCA_Total_score,PATNO,Elapsed_years,BL_Lymphocytes,Age,LEDD_sum,Disease_duration_diagnosis,Sex,base_UPSIT_total,base_RBDSQ_total_score)

      IDS <- c("PATNO","Elapsed_years")
      varIndex <- colnames(imp)

      impute.out <- amelia(as.data.frame(imp[, varIndex]), idvars = IDS,noms = c("Sex"),m = 100)
names(impute.out$imputations)<-c(1:100)

      Partner<-c()
      c<-c()
      imp_num<-c()
      interaction_effect_table_BL_Lymphocytes_sensitivity_imputation<-c()
      for(i in seq(1,100)){
        a1<-impute.out$imputations[i] %>% as.data.frame();colnames(a1)<-colnames(imp);
        a1 %>% dplyr::select(MOCA_Total_score,PATNO,Elapsed_years,BL_Lymphocytes,Age,LEDD_sum,Disease_duration_diagnosis,Sex,base_UPSIT_total,base_RBDSQ_total_score)->a1;
        a1$MOCA_Total_score<-scale(a1$MOCA_Total_score,center=T,scale=T);
        a1$BL_Lymphocytes<-scale(a1$BL_Lymphocytes,center=T,scale=T);
        a1$LEDD_sum<-scale(a1$LEDD_sum,center=T,scale=T);
        a1$Elapsed_years<-scale(a1$Elapsed_years,center=T,scale=T);
        a1$Disease_duration_diagnosis<-scale(a1$Disease_duration_diagnosis,center=T,scale=T);
        a1$base_UPSIT_total<-scale(a1$base_UPSIT_total,center=T,scale=T);
        a1$base_RBDSQ_total_score<-scale(a1$base_RBDSQ_total_score,center=T,scale=T);
        lme4::lmer(MOCA_Total_score~Elapsed_years*BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b
        lme4::lmer(MOCA_Total_score~Elapsed_years+BL_Lymphocytes+Age+LEDD_sum+Disease_duration_diagnosis+Sex+base_UPSIT_total+base_RBDSQ_total_score+(1|PATNO),a1,REML=FALSE)->b.null
        summary(b)->b.summary;confint(b) %>% suppressMessages() ->b.confint;
        b.summary$coefficients[10,]->e;b.confint[12,]->f;b.summary$AICtab[1:2]->A.BIC;
        anova(b,b.null)->d;summary(b.null)->d.summary;d.summary$AICtab[1:2]-b.summary$AICtab[1:2]->A.BIC
        c(e[1],e[2],e[3],f[1],f[2],d$`Pr(>Chisq)`[2],A.BIC)->d;
        bind_rows(c,d) %>% suppressMessages()->c;
        append(Partner,i)->Partner
      }
      colnames(c)[6]<-"p.value"
      interaction_effect_table_BL_Lymphocytes_sensitivity_imputation<-c


      rm(a1,b,b.confint,b.null,b.summary,c,d.summary,imp,impute.out,A.BIC,d,e,f,i,IDS,imp_num,Partner,varIndex)

