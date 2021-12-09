#ifndef TLATTICE_H
#define TLATTICE_H
#include <memory>
#include "TParticle.h"
#include <vector>
#include "RandUtils.h"
#include <SFML/Graphics.hpp>
#include "matrix.h"
#include "small_uint_set.hpp"

//"class forward declaration" da mettere quando 2 classi si includono a vicenda
class TParticle;

//! This class ...
class TLattice : public sf::Drawable {

public:
  //! Constructor
  TLattice();
  virtual ~TLattice();

  //Class Members
  //! Total number of particles stored in the Parts vector
  int N;
  //! Number of initial free particles
  int Nfree;
  //! Number of particles in the polimer
  int Nfix;
  //! Number of active free particles
  int Na;
  //! Number of passive free particles
  int Np;
  //! Fraction of activated particle: when overcomes 1, one particle is activated.
  double Activating;
  //! Max number of particles in the aggregate
  inline static int MAX_Nfix;
  //! Threshold for central sites activation
  inline static double act;
  //! True if a particles of the aggregates is on the border
  bool OutofGrid;
  //! The Lattice own as a member the vector of All Parts
  std::vector<TParticle> Parts;
  //! Number of YL and DL links
  int nYL, nDL;

  //Class Methods
  //! Fill the Lattice whit pN Particles
  void RandomFill(int pN);
    //! Set for a MF simulation
    void SetForMF();
  //! Set for a DLA simulation
  void SetForDLA();

  //! Make a time step forwar for all the system
  bool Evolve();
  bool EvolveMF();
  bool EvolveMF2();
  //! The Lattice is composed by Grid Elements made up with an int (the Index of the Particle if present) and the bool Is_Central
  //typedef mtx::small_vector<uint16_t> GridElement; // type def
  typedef mtx::small_uint_set<uint16_t,100> GridElement;

  //Class Setter ang Getter Methods
  //! Put pIndex of the Index-Particle in the pSite.x and pSite.y position of the Lattice
  void SetSitePosition(TSite &pSite, int pIndex);
  //! Clear the pSite.x and pSite.y position of the Lattice (at the moment, pIndex is unused)
  void ClearSitePosition(TSite &pSite, int pIndex);
  //! Get the vector of indexes contained in the position of pSite
  GridElement &GetSiteIndexes(TSite &pSite);
  //! Get the particle with a site in pSite (if there is any particle here)
  TParticle &GetParticle(int pIndex);

protected:
  void draw(sf::RenderTarget& target, sf::RenderStates states) const;

private:

  //! Matrix class manage positions of particles trought their indexes
  // std::array<std::array<GridElement, TSite::Ly>, TSite::Lx> Grid;
   mtx::matrix<GridElement, TSite::Lx, TSite::Ly> Grid;

  //! Add a particle in a specific position and spin
    void PutParticle(TSite pSite, int pSpin);

    //! Possibility of X coordinates translation on a triangular lattice
    static constexpr int dx[6] = {2, 1, -1, -2, -1, 1};  //IN ORDINE!!
    //! Possibility of Y coordinates translation on a triangular lattice
    static constexpr int dy[6] = {0, 1, 1, 0, -1, -1};

    //! Calculate probability to get a free monomer with active central sites
    double probAct();

    //! Inhibition (due to ratio between real and model interaction radius (10) or areas (100) )
    int Inh = 10;

};

#endif // TLATTICE_H
