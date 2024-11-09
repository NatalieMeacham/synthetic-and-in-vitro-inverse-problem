%S500
concvecSs500=[1     2     3     4     5     6     7     8     9    10    11];
evvecs500=[0.391554778598877
   0.352425300405143
   0.389966196518857
   0.474806025021250
   0.653831925103922
   0.827238854301867
   0.877074411165618
   1.000000000000000
   0.992086312797904
   1.000000000000000
   0.910351624119585];
concvecSs100=[ 1     2     3     4     5     6     7     8     9    10];
evvecs100=[0.309137258143782
   0.376226119126448
   0.444147924964657
   0.488238161747171
   0.575640508475282
   0.692432499975030
   0.820227387509498
   1.000000000000000
   1.000000000000000
   0.982673013036324];
concvecSr250=[1     2     3     4     5     6     7     8     9    10    11];
evvecr250=[0.002471250675863
   0.127051505876649
   0.134933309964287
   0.135628018300982
   0.183493012158774
   0.211724440249231
   0.149751623117030
   0.114473971160524
   0.000000000000000
   0.000000000000000
   0.023368902402449];
concvecSr500=[ 1     2     3     4     5     6     7     8     9    10    11];
evvecr500=[0.005653121488685
   0.030778546248698
   0.028191744271588
   0.000671181988044
   0.027802537570680
   0.046127126386338
   0.000000000000000
   0.000000000000000
   0.000000000000000
   0.000000000000000
   0.015621836156726];
concvecSbf11=[  1     2     3     4     5     6     7     8     9    10    11];
evvecbf11=[0.199500419290343
   0.270224097312327
   0.299152332415836
   0.295285747402492
   0.319495739999136
   0.367514320329439
   0.384949307641148
   0.369222334192770
   0.251851433123408
   0.170497317171638
   0.307429525772766];
concvecSbf12=[ 1     2     3     4     5     6     7     8     9    10    11];
evvecbf12=[0.122159159646169
   0.109959944033354
   0.121343177525710
   0.138427504349731
   0.184550933157761
   0.180404249604588
   0.227362252986868
   0.128476985398751
   0.069412681263310
   0.000000000000000
   0.138038080485024];
concvecSbf21=[1     2     3     4     5     6     7     8     9    10    11
];
evvecbf21=[0.392641562612595
   0.414362693626728
   0.424092262316530
   0.542558491547363
   0.626085478998792
   0.676012738599930
   0.727062601862941
   0.722953321626847
   0.653787072273230
   0.605415698766361
   0.696717225937701];
concvecSbf41=[1     2     3     4     5     6     7     8     9    10    11
];
evvecbf41=[ 0.338227965860502
   0.393182636770071
   0.413165524751177
   0.484858984341095
   0.618003179229237
   0.657787331045966
   0.769264286272438
   0.816442865706253
   0.737992224323920
   0.743279011796481
   0.875163598277047];

clr=jet(6);

dosevec= [0, 0.0260,0.0536, 0.1094, 0.2244, 0.3410,0.4603,1.1908,2.4430,3.7133,5.0119];
% 
% figure
% plot(concvecSs100,evvecs100,'d','LineWidth',2,'Color',clr(3,:))
% hold on
% plot(concvecSs500,evvecs500,'d','LineWidth',2,'Color',clr(3,:))

%just resistant
figure
semilogx(dosevec(concvecSr250),evvecr250,'d','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))
hold on
semilogx(dosevec(concvecSr500),evvecr500,'o','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))
% hold on
% semilogx(dosevec(concvecSs500),evvecs500,'d','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
% hold on
% semilogx(dosevec(concvecSs100),evvecs100,'d','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
% %legend('Resistant Data 1','Resistant Data 2','Sensitive Data 1','Sensitive Data 2','Location','northwest')
% semilogx(dosevec(concvecSbf12),evvecbf12,'d','MarkerSize',8,'LineWidth',2,'Color',clr(2,:))
% hold on
% semilogx(dosevec(concvecSbf11),evvecbf11,'d','MarkerSize',8,'LineWidth',2,'Color',clr(3,:))
% hold on
% semilogx(dosevec(concvecSbf21),evvecbf21,'d','MarkerSize',8,'LineWidth',2,'Color',clr(4,:))
% hold on
% semilogx(dosevec(concvecSbf41),evvecbf41,'d','MarkerSize',8,'LineWidth',2,'Color',clr(5,:))
legend('Resistant 1','Resistant 2','Location','northwest','FontSize',16)
ylim([0 1])
xlabel('Log(Dosage)')
ylabel('Mean Recovered Sensitivity')
set(gca,"FontSize",20)

%resistant and sensitive
figure
semilogx(dosevec(concvecSr250),evvecr250,'d','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))
hold on
semilogx(dosevec(concvecSr500),evvecr500,'o','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))
hold on
semilogx(dosevec(concvecSs500),evvecs500,'d','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
hold on
semilogx(dosevec(concvecSs100),evvecs100,'o','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
% %legend('Resistant Data 1','Resistant Data 2','Sensitive Data 1','Sensitive Data 2','Location','northwest')
% semilogx(dosevec(concvecSbf12),evvecbf12,'d','MarkerSize',8,'LineWidth',2,'Color',clr(2,:))
% hold on
% semilogx(dosevec(concvecSbf11),evvecbf11,'d','MarkerSize',8,'LineWidth',2,'Color',clr(3,:))
% hold on
% semilogx(dosevec(concvecSbf21),evvecbf21,'d','MarkerSize',8,'LineWidth',2,'Color',clr(4,:))
% hold on
% semilogx(dosevec(concvecSbf41),evvecbf41,'d','MarkerSize',8,'LineWidth',2,'Color',clr(5,:))
legend('Resistant 1','Resistant 2','Sensitive 1','Sensitive 2','Location','northwest','FontSize',16)
xlabel('Log(Dosage)')
ylabel('Mean Recovered Sensitivity')
set(gca,"FontSize",20)

%resistant, sensitive, and 1:2
figure
semilogx(dosevec(concvecSr250),evvecr250,'d','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))
hold on
semilogx(dosevec(concvecSr500),evvecr500,'o','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))
%legend('Resistant Data 1','Resistant Data 2','Sensitive Data 1','Sensitive Data 2','Location','northwest')
semilogx(dosevec(concvecSbf12),evvecbf12,'d','MarkerSize',8,'LineWidth',2,'Color',clr(2,:))
% hold on
% semilogx(dosevec(concvecSbf11),evvecbf11,'d','MarkerSize',8,'LineWidth',2,'Color',clr(3,:))
% hold on
% semilogx(dosevec(concvecSbf21),evvecbf21,'d','MarkerSize',8,'LineWidth',2,'Color',clr(4,:))
% hold on
% semilogx(dosevec(concvecSbf41),evvecbf41,'d','MarkerSize',8,'LineWidth',2,'Color',clr(5,:))
% hold on
semilogx(dosevec(concvecSs500),evvecs500,'d','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
hold on
semilogx(dosevec(concvecSs100),evvecs100,'o','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
legend('Resistant 1','Resistant 2','S:2R Mixture','Sensitive 1','Sensitive 2','Location','northwest','FontSize',16)
xlabel('Log(Dosage)')
ylabel('Mean Recovered Sensitivity')
set(gca,"FontSize",20)

%r,s,1:2,1:1
figure
semilogx(dosevec(concvecSr250),evvecr250,'d','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))
hold on
semilogx(dosevec(concvecSr500),evvecr500,'o','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))

%legend('Resistant Data 1','Resistant Data 2','Sensitive Data 1','Sensitive Data 2','Location','northwest')
semilogx(dosevec(concvecSbf12),evvecbf12,'d','MarkerSize',8,'LineWidth',2,'Color',clr(2,:))
hold on
semilogx(dosevec(concvecSbf11),evvecbf11,'d','MarkerSize',8,'LineWidth',2,'Color',clr(3,:))
% hold on
% semilogx(dosevec(concvecSbf21),evvecbf21,'d','MarkerSize',8,'LineWidth',2,'Color',clr(4,:))
% hold on
% semilogx(dosevec(concvecSbf41),evvecbf41,'d','MarkerSize',8,'LineWidth',2,'Color',clr(5,:))
% hold on
semilogx(dosevec(concvecSs500),evvecs500,'d','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
hold on
semilogx(dosevec(concvecSs100),evvecs100,'o','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
legend('Resistant 1','Resistant 2','S:2R Mixture','S:R Mixture','Sensitive 1','Sensitive 2','Location','northwest','FontSize',16)
xlabel('Log(Dosage)')
ylabel('Mean Recovered Sensitivity')
set(gca,"FontSize",20)

%all but 4:1
figure
semilogx(dosevec(concvecSr250),evvecr250,'d','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))
hold on
semilogx(dosevec(concvecSr500),evvecr500,'o','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))

%legend('Resistant Data 1','Resistant Data 2','Sensitive Data 1','Sensitive Data 2','Location','northwest')
semilogx(dosevec(concvecSbf12),evvecbf12,'d','MarkerSize',8,'LineWidth',2,'Color',clr(2,:))
hold on
semilogx(dosevec(concvecSbf11),evvecbf11,'d','MarkerSize',8,'LineWidth',2,'Color',clr(3,:))
hold on
semilogx(dosevec(concvecSbf21),evvecbf21,'d','MarkerSize',8,'LineWidth',2,'Color',clr(4,:))
% hold on
% semilogx(dosevec(concvecSbf41),evvecbf41,'d','MarkerSize',8,'LineWidth',2,'Color',clr(5,:))
hold on
semilogx(dosevec(concvecSs500),evvecs500,'d','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
hold on
semilogx(dosevec(concvecSs100),evvecs100,'o','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
legend('Resistant 1','Resistant 2','S:2R Mixture','S:R Mixture','2S:R Mixture','Sensitive 1','Sensitive 2','Location','northwest','FontSize',16)
xlabel('Log(Dosage)')
ylabel('Mean Recovered Sensitivity')
set(gca,"FontSize",20)


%all
figure
semilogx(dosevec(concvecSr250),evvecr250,'d','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))
hold on
semilogx(dosevec(concvecSr500),evvecr500,'o','MarkerSize',8,'LineWidth',2,'Color',clr(1,:))

%legend('Resistant Data 1','Resistant Data 2','Sensitive Data 1','Sensitive Data 2','Location','northwest')
semilogx(dosevec(concvecSbf12),evvecbf12,'d','MarkerSize',8,'LineWidth',2,'Color',clr(2,:))
hold on
semilogx(dosevec(concvecSbf11),evvecbf11,'d','MarkerSize',8,'LineWidth',2,'Color',clr(3,:))
hold on
semilogx(dosevec(concvecSbf21),evvecbf21,'d','MarkerSize',8,'LineWidth',2,'Color',clr(4,:))
hold on
semilogx(dosevec(concvecSbf41),evvecbf41,'d','MarkerSize',8,'LineWidth',2,'Color',clr(5,:))
hold on
semilogx(dosevec(concvecSs500),evvecs500,'d','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
hold on
semilogx(dosevec(concvecSs100),evvecs100,'o','MarkerSize',8,'LineWidth',2,'Color',clr(6,:))
legend('Resistant 1','Resistant 2','S:2R Mixture','S:R Mixture','2S:R Mixture','4S:R Mixture','Sensitive 1','Sensitive 2','Location','northwest','FontSize',16)
xlabel('Log(Dosage)')
ylabel('Mean Recovered Sensitivity')
set(gca,"FontSize",20)

