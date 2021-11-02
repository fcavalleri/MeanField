#include "TLattice.h"

static constexpr int Lx = TSite::Lx, Ly = TSite::Ly;

TLattice::TLattice() :  /*Grid(Lx, std::vector<GridElement>(Ly)),*/
N(0), Nfix(0), Nfree(0), Na(0), Np(0), nYL(0), nDL(0), OutofGrid(false), Activating(0) //Initialize
{
  // Set static parameters for TParticle and TSite classes
  TParticle::Lattice = this; //this è il puntatore all'istanza corrente della classe
}

TLattice::~TLattice() {
  //dtor
}

void TLattice::RandomFill(int pN) {

//Fill with a cycle the vector Parts
  for (int i = 0; i < pN; i++) {
    TParticle x{N};        //Create a Tparticle-Type variable named x Setting x Index as current N
    N++;                //Increment N by unit
    x.SetParticlePosition();
    Parts.push_back(x); //Put x in the last position of vector Parts
  }
}

void TLattice::PutParticle(TSite pSite, int pSpin) {
    TParticle x(Nfix, pSite, pSpin);        //Create a Tparticle-Type variable named x Setting x Index as current N
    Nfix++;                //Increment N by unit
    x.SetParticlePosition();
    Parts.push_back(x); //Put x in the last position of vector Parts
}

void TLattice::SetForMF() {

    Np = Nfree - 1;
    Na = 1;

    PutParticle(TSite(Lx / 2, Ly / 2), 0 );

    Parts[0].is_freeL = true;
    Parts[0].is_freeR = true;
    Parts[0].is_activeA = true;
    Parts[0].is_activeB = true;
    Parts[0].mob = TParticle::MobState::BLOCKED;

    Parts[0].SetParticlePosition();
}

void TLattice::SetForDLA() {

  Parts[0].ClearParticlePosition();

  Parts[0].CSite.x = Lx / 2;
  Parts[0].CSite.y = Ly / 2;
  Parts[0].Spin = 0;
  Parts[0].LSite.x = Lx / 2 - 2;
  Parts[0].LSite.y = Ly / 2;
  Parts[0].RSite.x = Lx / 2 + 2;
  Parts[0].RSite.y = Ly / 2;
  Parts[0].is_freeL = true;
  Parts[0].is_freeR = true;
  Parts[0].is_activeA = true;
  Parts[0].is_activeB = true;
  Parts[0].mob = TParticle::MobState::BLOCKED;

  Parts[0].SetParticlePosition();

  Nfix=1;
}

bool TLattice::Evolve() {
  //For N times, chose a random particle from Parts vector and Evolve it
  for (int i = 0; i < N; i++) {
      int j=randM(N);
    if (Parts[j].Evolve()) {
        // Particles in polimer
        Nfix++;
        if (Nfix > (MAX_Nfix-1) || OutofGrid) return true;
        //Add new particle to preserve free monomer concentration
        //RandomFill(1);
    };
  }
    return false;
}


bool TLattice::EvolveMF2() {

    bool result=true;

    // Tengo gli indici dei nuovi monomeri creati:
    std::vector<int> IndexCreatedParticles;

    // Calcolo la probabilità di estrarre un monomero libero con i siti centrali attivi
    double pA = probAct();

    int NfixMemory = Nfix;  // per evitare che la condizione di ciclo fos si sposti aggiornando Nfix

        for (int i = 0; i < NfixMemory; i++) {

            int j;  //Indice particella su cui simulare un nuovo attacco o chiusura
            if (Nfix == 1) { j = 0; } else { j = randM(Nfix); }

            // Evolvo con il metodo già costruito: se BLOCKED non succede niente, se LINKED prova a chiudersi.
            if (Parts[j].Evolve())
                if (Nfix > (MAX_Nfix - 1) || OutofGrid)
                    return true;  // se si chiude controlla la griglia

            // Provo a vedere se ci sono le condizioni per creare ed attaccare un nuovo monomero:

            if (ranMT() < pA) { // nuovo possibile monomero Attivo (siti centrali già attivati)

                // Ci sono 6 modi diversi in cui si può attaccare un nuovo monomero con i siti A e B attivi: li testo tutti

                // DLAs
                if (Parts[j].LinkedWith[2] == -1 && Parts[j].LinkedWith[3] == -1 && Parts[j].is_activeA) {

                    if (ranMT() < 2 * (double) (Nfree - Nfix) / (6 * Lx * Ly)) {
                        // Add a particle in the polimer
                        PutParticle(Parts[j].RSite, (Parts[j].Spin + 3) % 6);
                        IndexCreatedParticles.push_back(Nfix - 1);
                        Parts[Nfix - 1].is_activeA = true;
                        Parts[Nfix - 1].is_activeB = true;
                        if (!Parts[Nfix - 1].Evolve()) throw std::runtime_error("Error:");
                    }
                }

                // YLA
                if (Parts[j].LinkedWith[2] == -1 && Parts[j].is_activeA) {

                    if (ranMT() < 2 * (double) (Nfree - Nfix) / (6 * Lx * Ly)) {
                        // Add a particle in the polimer
                        PutParticle(
                                Parts[j].CSite.GetTranslatedCSite(dx[(Parts[j].Spin + 1) % 6],
                                                                  dy[(Parts[j].Spin + 1) % 6]),
                                (Parts[j].Spin + 4) % 6);
                        IndexCreatedParticles.push_back(Nfix - 1);
                        Parts[Nfix - 1].is_activeA = true;
                        Parts[Nfix - 1].is_activeB = true;
                        if (!Parts[Nfix - 1].Evolve()) throw std::runtime_error("Error:");
                    }
                }

                // DLBs
                if (Parts[j].LinkedWith[0] == -1 && Parts[j].LinkedWith[1] == -1 && Parts[j].is_activeB) {

                    if (ranMT() < 2 * (double) (Nfree - Nfix) / (6 * Lx * Ly)) {
                        // Add a particle in the polimer
                        PutParticle(Parts[j].LSite, (Parts[j].Spin + 3) % 6);
                        IndexCreatedParticles.push_back(Nfix - 1);
                        Parts[Nfix - 1].is_activeA = true;
                        Parts[Nfix - 1].is_activeB = true;
                        if (!Parts[Nfix - 1].Evolve()) throw std::runtime_error("Error:");
                    }
                };


                // YLB
                if (Parts[j].LinkedWith[1] == -1 && Parts[j].is_activeB) {

                    if (ranMT() < 2 * (double) (Nfree - Nfix) / (6 * Lx * Ly)) {
                        // Add a particle in the polimer
                        PutParticle(
                                Parts[j].CSite.GetTranslatedCSite(dx[(Parts[j].Spin + 2) % 6],
                                                                  dy[(Parts[j].Spin + 2) % 6]),
                                (Parts[j].Spin + 2) % 6);
                        IndexCreatedParticles.push_back(Nfix - 1);
                        Parts[Nfix - 1].is_activeA = true;
                        Parts[Nfix - 1].is_activeB = true;
                        if (!Parts[Nfix - 1].Evolve()) throw std::runtime_error("Error:");
                    }
                };

                // YLR
                if (Parts[j].LinkedWith[3] == -1) {

                    if (ranMT() < 2 * (double) (Nfree - Nfix) / (6 * Lx * Ly)) {
                        // Add a particle in the polimer
                        PutParticle(Parts[j].RSite, (Parts[j].Spin + 2) % 6);
                        IndexCreatedParticles.push_back(Nfix - 1);
                        Parts[Nfix - 1].is_activeA = true;
                        Parts[Nfix - 1].is_activeB = true;
                        if (!Parts[Nfix - 1].Evolve()) throw std::runtime_error("Error:");
                    }
                };

                // YLL
                if (Parts[j].LinkedWith[0] == -1) {

                    if (ranMT() < 2 * (double) (Nfree - Nfix) / (6 * Lx * Ly)) {
                        // Add a particle in the polimer
                        PutParticle(Parts[j].LSite, (Parts[j].Spin + 4) % 6);
                        IndexCreatedParticles.push_back(Nfix - 1);
                        Parts[Nfix - 1].is_activeA = true;
                        Parts[Nfix - 1].is_activeB = true;
                        if (!Parts[Nfix - 1].Evolve()) throw std::runtime_error("Error:");
                    }
                };

            } else { // nuovo possibile monomero passivo (siti centrali non ancora attivati)

                // Ci sono solo 2 modi diversi in cui si può attaccare un nuovo monomero passivo:

                // YLA
                if (Parts[j].LinkedWith[2] == -1 && Parts[j].is_activeA ) {

                    if (ranMT() < 2 * (double) (Nfree - Nfix) / (6 * Lx * Ly)) {
                        // Add a particle in the polimer
                        PutParticle(
                                Parts[j].CSite.GetTranslatedCSite(dx[(Parts[j].Spin + 1) % 6],
                                                                  dy[(Parts[j].Spin + 1) % 6]),
                                (Parts[j].Spin + 4) % 6);
                        IndexCreatedParticles.push_back(Nfix - 1);
                        if (!Parts[Nfix - 1].Evolve()) throw std::runtime_error("Error:");
                    }
                }

                // YLB
                if (Parts[j].LinkedWith[1] == -1 && Parts[j].is_activeB) {

                    if (ranMT() < 2 * (double) (Nfree - Nfix) / (6 * Lx * Ly)) {
                        // Add a particle in the polimer
                        PutParticle(
                                Parts[j].CSite.GetTranslatedCSite(dx[(Parts[j].Spin + 2) % 6],
                                                                  dy[(Parts[j].Spin + 2) % 6]),
                                (Parts[j].Spin + 2) % 6);
                        IndexCreatedParticles.push_back(Nfix - 1);
                        if (!Parts[Nfix - 1].Evolve()) throw std::runtime_error("Error:");
                    }
                }

            }
        }

    if (Nfix > (MAX_Nfix-1) || OutofGrid) return true;

    // Check if all created particles where correctly linked in the polimer
    for (int i : IndexCreatedParticles)
        if (Parts[i].mob == TParticle::MobState::FREE)
           std::cout << "WRONGCONDITION!\n";
            // arrivano sempre particelle con i siti centrali attivati!! qualcosa non va.
            // ne arrivano anche senza passare dal caso di particelle inattive (gli altri breakpoints)
        return false;
}

bool TLattice::EvolveMF() {
    //For N times, chose a random particle from Parts vector and Evolve it
    for (int i = 0; i < Nfix; i++) {

        int j;
        if (Nfix==1) {j=0;} else {j=randM(Nfix);}

        // Prova a chiudersi:
       if (Parts[j].mob==TParticle::MobState::LINKED) {

            if (Parts[j].CheckClose()) {
                Parts[j].CheckBorder();
                break;
            }
        }

 /*       // DLAs
        if (Parts[j].LinkedWith[2]==-1 && Parts[j].LinkedWith[3]==-1 ) {

            if (ranMT()<2*(double)(Nfree-Nfix)/(6*Lx*Ly)) {
                // Add a particle in the polimer
                PutParticle(Parts[j].RSite,(Parts[j].Spin+3)%6);
                Parts[Nfix-1].DLAs(Parts[j]);
                Parts[Nfix-1].CheckBorder();
            }
        }

        // YLA
        if (Parts[j].LinkedWith[2]==-1) {

            if (ranMT()<2*(double)(Nfree-Nfix)/(6*Lx*Ly)) {
                // Add a particle in the polimer
                PutParticle(Parts[j].CSite.GetTranslatedCSite(dx[(Parts[j].Spin + 1) % 6], dy[(Parts[j].Spin + 1) % 6]),(Parts[j].Spin+4)%6);
                Parts[Nfix-1].YLA(Parts[j]);
                Parts[Nfix-1].CheckBorder();
            }
        }

        // DLBs
        if (Parts[j].LinkedWith[0]==-1 && Parts[j].LinkedWith[1]==-1 ) {

            if (ranMT()<2*(double)(Nfree-Nfix)/(6*Lx*Ly)) {
                // Add a particle in the polimer
                PutParticle(Parts[j].LSite,(Parts[j].Spin+3)%6);
                Parts[Nfix-1].DLBs(Parts[j]);
                Parts[Nfix-1].CheckBorder();
            }
        };


        // YLB
        if (Parts[j].LinkedWith[1]==-1) {

            if (ranMT()<2*(double)(Nfree-Nfix)/(6*Lx*Ly)) {
                // Add a particle in the polimer
                PutParticle(Parts[j].CSite.GetTranslatedCSite(dx[(Parts[j].Spin + 2) % 6], dy[(Parts[j].Spin + 2) % 6]),(Parts[j].Spin+2)%6);
                Parts[Nfix-1].YLB(Parts[j]);
                Parts[Nfix-1].CheckBorder();
            }
        };
*/

        // YLR
        if (Parts[j].LinkedWith[3]==-1) {

            if (ranMT()<2*(double)(Nfree-Nfix)/(6*Lx*Ly)) {
                // Add a particle in the polimer
                PutParticle(Parts[j].RSite,(Parts[j].Spin+2)%6);
                Parts[Nfix-1].YLR(Parts[j]);
                Parts[Nfix-1].CheckBorder();
            }
        };

        // YLL
        if (Parts[j].LinkedWith[0]==-1) {

            if (ranMT()<2*(double)(Nfree-Nfix)/(6*Lx*Ly)) {
                // Add a particle in the polimer
                PutParticle(Parts[j].LSite,(Parts[j].Spin+4)%6);
                Parts[Nfix-1].YLL(Parts[j]);
                Parts[Nfix-1].CheckBorder();
            }
        };

        if (Nfix > (MAX_Nfix-1) || OutofGrid) return true;
    }
    return false;
}

void TLattice::SetSitePosition(TSite &pSite, int pIndex) {
  //Gives at the int element of the Grid struct the value pIndex
  Grid[pSite.x][pSite.y].insert_safe(pIndex);
  //Gives at the bool element of the Grid struct the value is_central given by pSite
  //Grid[pSite.x*L+pSite.y].is_central=pSite.is_central;
}

void TLattice::ClearSitePosition(TSite &pSite, int pIndex) {
  //Gives at the int element of the Grid struct the value -1. The bool value doesn't matter
  GridElement &cell = Grid[pSite.x][pSite.y];
  cell.erase(pIndex);
  //prima era
  /*std::remove(Grid[pos].begin(),Grid[pos].end(),pIndex),Grid[pos].end()*/
}

TLattice::GridElement &TLattice::GetSiteIndexes(TSite &pSite) {
//if (Grid[pSite.x*L+pSite.y].size()==0) return -1
  return Grid[pSite.x][pSite.y];
}

TParticle &TLattice::GetParticle(int pIndex) {
//not a copy, the reference at that particle with index pIndex
  return Parts.at(pIndex);
}

void TLattice::draw(sf::RenderTarget &target, sf::RenderStates states) const {
  for (auto &i : Parts) {
    if (i.mob != TParticle::MobState::FREE) {
    //if (!(i.CSite.x < 2 || i.CSite.x > (Lx- 3) || i.CSite.y < 1 || i.CSite.y > (Ly-2) )) {
          sf::Vertex monomer[] = {sf::Vertex(sf::Vector2f(i.LSite.x, i.LSite.y)),
                                  sf::Vertex(sf::Vector2f(i.RSite.x, i.RSite.y))};
          target.draw(monomer, 2, sf::Lines, states);
      }
    // Uncommento to draw also free monomers
    else {
        sf::Vertex monomerCM=sf::Vertex(sf::Vector2f(i.CSite.x,i.CSite.y), sf::Color::White);
        target.draw(&monomerCM,1,sf::Points);
    }
    //}
  }
}

double TLattice::probAct() {

    double pA;

    // Incremento di un fattore proporzionale al numero di monomeri inattivi * la probabilità di attivazione
    Activating += Np * act;

    // Aumento il di particelle attive quando la frazione raggiunge un intero
    int IntActivating= floor(Activating);
    Na = Na + IntActivating;
    Np = Np - IntActivating;

    // Accumulo in fattore restante per i successivi incrementi
    Activating -= IntActivating;

    // Probabilità di estrarre un monomero libero già attivo

    return pA = double(Na) / (Nfree);

    // TODO: conteggiare solo le particelle FREE
}

