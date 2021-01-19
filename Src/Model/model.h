// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QSet>
#include <QTextStream>

namespace LoboLab {

class ModelProd;
class ModelLink;

class Model {

 public:
  Model();
  ~Model();

  Model(const Model &source);
  Model &operator=(const Model &source);

  static Model *createRandom(int nProducts, int nMorphogens);
  static void cross(const Model *parent1, const Model *parent2,
                    Model *&child1, Model *&child2, int nMorphogens,
                    int nMinProducts, int nMaxProducts, int nMaxLinks);

  QSet<int> calcProductLabels() const;
  QSet<int> calcProductLabelsInUse(int nMorphogens) const;

  inline int nProducts() const { return products_.size(); }
  inline ModelProd *product(int i) const { return products_.at(i); }
  ModelProd *prodWithLabel(int label, int *i = NULL) const;
  
  int addRandomProduct();
  void addRandomProduct(int label, int nMorphogens);
  void duplicateProduct(int i, int nMorphogens);
  void removeProduct(int i);
  void removeProductWithLabel(int label);

  inline int nLinks() const { return links_.size(); }
  inline ModelLink *link(int i) const { return links_.at(i); }
  inline QList<ModelLink*> links() const { return links_; }
  QList<ModelLink*> linksToLabel(int label) const;
  QList<ModelLink*> calcLinksInUse(int nStrucLabels) const;
  int calcNLinksFromProd(int label) const;
  ModelLink *findLink(int regulator, int regulated) const;

  void addOrReplaceRandomLink();
  void addOrReplaceRandomLink(int regulatorLabel, int regulatedLabel);
  void duplicateLink(int i, int nMorphogens);
  void removeLink(int i);
  void removeLink(ModelLink *modelLink);
  void removeLink(int regulatorLabel, int regulatedLabel);

  int calcComplexity() const;
  int calcComplexityInUse(int nStrucProd) const;

  void mutate(int nMorphogens, int nMinProducts, int nMaxProducts, int nMaxLinks);
  void removeExcess(int nMinProducts, int nMaxProducts, int nMaxLinks);
  void clear();
  
  // Text Serialization
  void loadFromString(QString &str);
  QString toString();

  friend QTextStream &operator<<(QTextStream &stream, const Model &model);
  friend QTextStream &operator>>(QTextStream &stream, Model &model);

  static double parseDouble(QTextStream &stream);

 private:
  static void distributeProducts(const QSet<int> &fromProds, 
                                 QSet<int> *toProds1, QSet<int> *toProds2);
  static void copyProductsAndLinks(Model *to,
                                const Model *from1, const QSet<int> &products1,
                                const Model *from2, const QSet<int> &products2,
								int numMorphogens);
  static void copyProducts(Model *to,const Model *from, 
                           const QSet<int> &products);
  static void copyLinks(Model *to, 
                        const Model *from1, const QSet<int> &products1,
                        const Model *from2, const QSet<int> &products2,
						int numMorphogens);
  QMap<int, QList<int> > calcProductRegulators() const;
  int createNewLabel();


  QList<ModelProd*> products_;
  QList<ModelLink*> links_;
};

} // namespace LoboLab
